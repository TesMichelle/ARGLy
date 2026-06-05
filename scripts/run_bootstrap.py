import sys
import os
import subprocess
import tskit
import msprime
import pandas as pd
import numpy as np
import shutil
from scipy.optimize import minimize

# Target demographic parameters to test (Scenario 2 estimates)
TRUE_T_MIG = 10000.0
TRUE_N_GHOST = 10000.0
FIXED_T_SEP = 30000.0
FIXED_T_INTRO = 50000.0

def simulate_locus(prefix, seed):
    """
    Simulates a 50 kb locus under the target parameters.
    """
    demography = msprime.Demography()
    demography.add_population(name="AFR", initial_size=10000)
    demography.add_population(name="SAM", initial_size=10000)
    demography.add_population(name="Ghost", initial_size=10000)
    demography.add_population(name="AFR_stay", initial_size=10000)
    demography.add_population(name="ANC", initial_size=10000)
    demography.add_population(name="ARCH", initial_size=10000)
    
    demography.add_population_split(time=FIXED_T_INTRO, derived=["Ghost", "ANC"], ancestral="ARCH")
    demography.add_population_split(time=FIXED_T_SEP, derived=["AFR_stay", "SAM"], ancestral="ANC")
    demography.add_admixture(time=TRUE_T_MIG, derived="AFR", ancestral=["AFR_stay", "Ghost"], proportions=[0.8, 0.2])
    demography.sort_events()
    
    ts = msprime.sim_ancestry(
        samples={"AFR": 10, "SAM": 10},
        demography=demography,
        sequence_length=50_000,
        recombination_rate=1e-8,
        random_seed=seed
    )
    mutated_ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=seed)
    
    vcf_path = f"{prefix}.vcf"
    with open(vcf_path, "w") as f:
        mutated_ts.write_vcf(f)
    return vcf_path

def prepare_for_mle(input_ts_path, output_ts_path):
    ts = tskit.load(input_ts_path)
    tables = ts.dump_tables()
    if len(tables.populations) == 0:
        tables.populations.add_row()
    num_nodes = len(tables.nodes)
    new_pops = tables.nodes.population.copy()
    for u in range(num_nodes):
        orig_pop = new_pops[u]
        if orig_pop == 1:  # SAM
            new_pops[u] = tskit.NULL
        else:
            new_pops[u] = 0
    tables.nodes.population = new_pops
    tables.sort()
    tables.tree_sequence().dump(output_ts_path)

def run_singer(vcf_path, output_ts_path, seed):
    vcf_prefix = os.path.splitext(vcf_path)[0]
    singer_bin = os.environ.get("SINGER_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/singer_native")
    cmd = [
        singer_bin, "-input", vcf_prefix, "-output", f"{vcf_prefix}_inferred",
        "-start", "0", "-end", "50000", "-n", "30", "-thin", "3",
        "-Ne", "10000", "-r", "1e-8", "-m", "1.2e-8", "-seed", str(seed)
    ]
    subprocess.run(cmd, capture_output=True, check=True)
    
    last_iteration = 27
    convert_cmd = [
        sys.executable,
        os.environ.get("CONVERT_TSKIT_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/convert_to_tskit"),
        "-input", f"{vcf_prefix}_inferred", "-output", f"{vcf_prefix}_inferred_converted",
        "-start", str(last_iteration), "-end", str(last_iteration + 1), "-step", "1"
    ]
    subprocess.run(convert_cmd, capture_output=True, check=True)
    shutil.copy(f"{vcf_prefix}_inferred_converted_{last_iteration}.trees", output_ts_path)

def run_official_polegon(vcf_prefix, output_ts_path, seed):
    singer_nodes = f"{vcf_prefix}_inferred_nodes_27.txt"
    singer_branches = f"{vcf_prefix}_inferred_branches_27.txt"
    singer_muts = f"{vcf_prefix}_inferred_muts_27.txt"
    
    shutil.copy(singer_nodes, f"{vcf_prefix}_inferred_nodes.txt")
    shutil.copy(singer_branches, f"{vcf_prefix}_inferred_branches.txt")
    shutil.copy(singer_muts, f"{vcf_prefix}_inferred_muts.txt")
    
    polegon_bin = os.environ.get("POLEGON_PATH", "/home/igor/ARG_Research/POLEGON/POLEGON/POLEGON/polegon_native")
    cmd = [
        polegon_bin, "-input", f"/home/igor/ARG_Research/{vcf_prefix}_inferred",
        "-Ne", "10000", "-m", "1.2e-8", "-burn_in", "5", "-num_samples", "10",
        "-thin", "2", "-seed", str(seed)
    ]
    subprocess.run(cmd, capture_output=True, check=True)
    
    shutil.copy(f"{vcf_prefix}_inferred_new_nodes.txt", f"{vcf_prefix}_inferred_nodes_999.txt")
    shutil.copy(singer_branches, f"{vcf_prefix}_inferred_branches_999.txt")
    shutil.copy(singer_muts, f"{vcf_prefix}_inferred_muts_999.txt")
    
    convert_cmd = [
        sys.executable,
        os.environ.get("CONVERT_TSKIT_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/convert_to_tskit"),
        "-input", f"{vcf_prefix}_inferred", "-output", f"{vcf_prefix}_inferred_optimized_converted",
        "-start", "999", "-end", "1000", "-step", "1"
    ]
    subprocess.run(convert_cmd, capture_output=True, check=True)
    shutil.copy(f"{vcf_prefix}_inferred_optimized_converted_999.trees", output_ts_path)

def evaluate_mle(ts_path, t_mig, n_ghost):
    arg_target = "my_arg(2).arg"
    if os.path.exists(arg_target):
        os.remove(arg_target)
    shutil.copy(ts_path, arg_target)
    
    mle_bin = "./builddir/main"
    if not os.path.exists(mle_bin):
        mle_bin = "../builddir/main"
        if not os.path.exists(mle_bin):
            mle_bin = "main"
    cmd = [mle_bin, f"{t_mig:.6f}", f"{n_ghost:.6f}", f"{FIXED_T_SEP:.6f}", f"{FIXED_T_INTRO:.6f}"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return float(res.stdout.strip())
    except Exception:
        return None

def optimize_replicate(ts_path, rep_idx):
    def loss_func(params):
        t_mig, n_ghost = params
        likelihood = evaluate_mle(ts_path, t_mig, n_ghost)
        return -likelihood if likelihood is not None else 1e10
        
    x0 = [12000.0, 12000.0]
    bounds = [(100.0, FIXED_T_SEP - 50.0), (1000.0, 50000.0)]
    res = minimize(loss_func, x0, method='Nelder-Mead', bounds=bounds, options={'maxiter': 35, 'xatol': 50.0, 'fatol': 5.0})
    return res.x

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, ".."))
    os.chdir(project_root)
    os.makedirs("temp_data", exist_ok=True)
    num_replicates = 100
    estimates = []
    
    print(f"=== LAUNCHING PARAMETRIC BOOTSTRAP ({num_replicates} Replicates) ===")
    
    for r in range(1, num_replicates + 1):
        seed = 99 + r * 19
        prefix = f"temp_data/boot_temp_rep{r}"
        
        try:
            # 1. Simulate 50 kb locus
            vcf_path = simulate_locus(prefix, seed)
            
            # 2. Run SINGER
            inferred_ts = f"{prefix}_inferred.trees"
            run_singer(vcf_path, inferred_ts, seed)
            
            # 3. Run POLEGON scaling
            optimized_ts = f"{prefix}_optimized.trees"
            run_official_polegon(prefix, optimized_ts, seed)
            
            # 4. Format population structure
            mle_ready_ts = f"{prefix}_mle_ready.trees"
            prepare_for_mle(optimized_ts, mle_ready_ts)
            
            # 5. Optimize demographic parameters for this replicate
            est_t_mig, est_n_ghost = optimize_replicate(mle_ready_ts, r)
            estimates.append({"Replicate": r, "Est_T_MIG": est_t_mig, "Est_N_GHOST": est_n_ghost})
            print(f"[Rep {r}/{num_replicates}] Completed: Est T_MIG={est_t_mig:.2f}, Est N_GHOST={est_n_ghost:.2f}")
            
        except Exception as e:
            print(f"[Rep {r}/{num_replicates}] Failed with error: {e}")
            
        finally:
            # Cleanup temporary files for this replicate
            for ext in [".vcf", "_inferred.trees", "_optimized.trees", "_mle_ready.trees", "_inferred_nodes.txt", "_inferred_branches.txt", "_inferred_muts.txt"]:
                f = f"{prefix}{ext}"
                if os.path.exists(f):
                    os.remove(f)
            # Remove converted files
            f_conv = f"{prefix}_inferred_converted_27.trees"
            if os.path.exists(f_conv): os.remove(f_conv)
            f_opt_conv = f"{prefix}_inferred_optimized_converted_999.trees"
            if os.path.exists(f_opt_conv): os.remove(f_opt_conv)
            
    df = pd.DataFrame(estimates)
    df.to_csv("bootstrap_estimates.csv", index=False)
    
    # Calculate summary statistics
    mean_t_mig = df["Est_T_MIG"].mean()
    se_t_mig = df["Est_T_MIG"].std()
    ci_lower_t_mig = df["Est_T_MIG"].quantile(0.025)
    ci_upper_t_mig = df["Est_T_MIG"].quantile(0.975)
    
    mean_n_ghost = df["Est_N_GHOST"].mean()
    se_n_ghost = df["Est_N_GHOST"].std()
    ci_lower_n_ghost = df["Est_N_GHOST"].quantile(0.025)
    ci_upper_n_ghost = df["Est_N_GHOST"].quantile(0.975)
    
    print("\n" + "="*80)
    print(f"{'BOOTSTRAP SUMMARY STATISTICS (B = 100)':^80}")
    print("="*80)
    print(f"T_MIG:   Mean = {mean_t_mig:.2f}, SE = {se_t_mig:.2f}")
    print(f"         95% CI (Percentile) = [{ci_lower_t_mig:.2f}, {ci_upper_t_mig:.2f}]")
    print("-"*80)
    print(f"N_GHOST: Mean = {mean_n_ghost:.2f}, SE = {se_n_ghost:.2f}")
    print(f"         95% CI (Percentile) = [{ci_lower_n_ghost:.2f}, {ci_upper_n_ghost:.2f}]")
    print("="*80)

if __name__ == "__main__":
    main()
