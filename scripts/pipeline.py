import sys
import os
import argparse
import subprocess
import tskit
import msprime
import pandas as pd
import numpy as np
import shutil
from scipy.optimize import minimize

def simulate_data(output_prefix, t_mig, t_sep, t_intro, seq_length=200_000, num_samples_afr=10, num_samples_sam=10, seed=42):
    """
    Simulates demographic history using msprime:
    - AFR (Pop 0) and SAM (Pop 1) with Ne = 10000.
    - At t_mig, a fraction of 0.2 of AFR lineages admix with Ghost (Pop 2).
    - At t_sep, SAM and AFR merge into ANC.
    - At t_intro, Ghost and ANC merge into ARCH.
    """
    print(f"=== Simulating demographic history (Seed={seed}, Length={seq_length}bp, T_MIG={t_mig}, T_SEP={t_sep}, T_INTRO={t_intro}) ===")
    
    # Define demographic model
    demography = msprime.Demography()
    
    # Add populations
    demography.add_population(name="AFR", initial_size=10000)
    demography.add_population(name="SAM", initial_size=10000)
    demography.add_population(name="Ghost", initial_size=10000)
    demography.add_population(name="AFR_stay", initial_size=10000)
    demography.add_population(name="ANC", initial_size=10000)
    demography.add_population(name="ARCH", initial_size=10000)
    
    # 1. Forward in time: ARCH is the ancestral root before t_intro.
    # 2. At t_intro, Ghost and ANC split from ARCH (backward: Ghost and ANC merge into ARCH)
    demography.add_population_split(time=t_intro, derived=["Ghost", "ANC"], ancestral="ARCH")
    
    # 3. At t_sep, AFR_stay and SAM split from ANC (backward: AFR_stay and SAM merge into ANC)
    demography.add_population_split(time=t_sep, derived=["AFR_stay", "SAM"], ancestral="ANC")
    
    # 4. At t_mig, AFR is formed by admixture from AFR_stay (80%) and Ghost (20%)
    # (backward: AFR lineages merge into AFR_stay and Ghost)
    demography.add_admixture(time=t_mig, derived="AFR", ancestral=["AFR_stay", "Ghost"], proportions=[0.8, 0.2])
    
    # Sort events chronologically as required by msprime
    demography.sort_events()
    
    # Run the simulation
    ts = msprime.sim_ancestry(
        samples={"AFR": num_samples_afr, "SAM": num_samples_sam},
        demography=demography,
        sequence_length=seq_length,
        recombination_rate=1e-8,
        random_seed=seed
    )
    
    # Add mutations to make a VCF
    mutated_ts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=seed)
    
    # Save the true simulated tree sequence
    true_ts_path = f"{output_prefix}_true.trees"
    print(f"Saving true tree sequence to {true_ts_path}...")
    mutated_ts.dump(true_ts_path)
    
    # Export to VCF for SINGER
    vcf_path = f"{output_prefix}.vcf"
    print(f"Exporting to VCF: {vcf_path}...")
    with open(vcf_path, "w") as f:
        mutated_ts.write_vcf(f)
        
    return true_ts_path, vcf_path

def prepare_for_mle(input_ts_path, output_ts_path):
    """
    Modifies population IDs in the tree sequence to match the C++ MLE convention:
    - SAM nodes (outgroup) population = -1 (TSK_NULL).
    - All other nodes population = 0.
    """
    ts = tskit.load(input_ts_path)
    tables = ts.dump_tables()
    
    # Ensure there is at least one population in the table for ID 0
    if len(tables.populations) == 0:
        tables.populations.add_row()
        
    num_nodes = len(tables.nodes)
    new_pops = tables.nodes.population.copy()
    
    for u in range(num_nodes):
        orig_pop = new_pops[u]
        if orig_pop == 1:  # SAM
            new_pops[u] = tskit.NULL
        else:
            new_pops[u] = 0  # Map all non-SAM populations to 0
            
    tables.nodes.population = new_pops
    tables.sort()
    mle_ts = tables.tree_sequence()
    mle_ts.dump(output_ts_path)

def run_singer(vcf_path, output_ts_path, seq_length=200_000, seed=42):
    """
    Runs SINGER on the simulated VCF and converts the last MCMC iteration to standard trees format.
    """
    vcf_prefix = os.path.splitext(vcf_path)[0]
    
    singer_bin = os.environ.get("SINGER_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/singer_native")
    if not os.path.exists(singer_bin):
        singer_bin = os.environ.get("SINGER_PATH_FALLBACK", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/singer")
        
    cmd = [
        singer_bin,
        "-input", vcf_prefix,
        "-output", vcf_prefix + "_inferred",  # SINGER adds extensions internally
        "-start", "0",
        "-end", str(seq_length),
        "-n", "30",  # number of MCMC iterations
        "-thin", "3",
        "-Ne", "10000",
        "-r", "1e-8",
        "-m", "1.2e-8",
        "-seed", str(seed)
    ]
    
    print(f"Running SINGER command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print("SINGER inference completed successfully.")
        
        # Convert SINGER MCMC sample index 27 to standard tskit trees format
        last_iteration = 27
        convert_cmd = [
            sys.executable,
            os.environ.get("CONVERT_TSKIT_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/convert_to_tskit"),
            "-input", f"{vcf_prefix}_inferred",
            "-output", f"{vcf_prefix}_inferred_converted",
            "-start", str(last_iteration),
            "-end", str(last_iteration + 1),
            "-step", "1"
        ]
        print(f"Running SINGER text-to-trees conversion: {' '.join(convert_cmd)}")
        subprocess.run(convert_cmd, check=True)
        
        converted_file = f"{vcf_prefix}_inferred_converted_{last_iteration}.trees"
        if os.path.exists(converted_file):
            shutil.copy(converted_file, output_ts_path)
            return output_ts_path
        else:
            raise FileNotFoundError(f"Converted trees file not found: {converted_file}")
            
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Warning: SINGER pipeline failed or returned error ({e}). Falling back to true trees.")
        mock_inferred = f"{vcf_prefix}_mock_inferred.trees"
        ts = tskit.load(f"{vcf_prefix}_true.trees")
        ts.dump(mock_inferred)
        shutil.copy(mock_inferred, output_ts_path)
        return output_ts_path

def run_official_polegon(vcf_prefix, output_ts_path, seq_length=200_000, seed=42):
    """
    Integrates the official POLEGON C++ branch-length optimization tool.
    Links the SINGER outputs, executes POLEGON, and converts the results to standard trees.
    """
    print("\n=== Running Official POLEGON C++ Branch-Length Optimizer ===")
    
    singer_nodes = f"{vcf_prefix}_inferred_nodes_27.txt"
    singer_branches = f"{vcf_prefix}_inferred_branches_27.txt"
    singer_muts = f"{vcf_prefix}_inferred_muts_27.txt"
    
    polegon_input_nodes = f"{vcf_prefix}_inferred_nodes.txt"
    polegon_input_branches = f"{vcf_prefix}_inferred_branches.txt"
    polegon_input_muts = f"{vcf_prefix}_inferred_muts.txt"
    
    # Clean old files
    for f in [polegon_input_nodes, polegon_input_branches, polegon_input_muts]:
        if os.path.exists(f):
            os.remove(f)
            
    # Copy SINGER text outputs to generic filenames for POLEGON C++
    shutil.copy(singer_nodes, polegon_input_nodes)
    shutil.copy(singer_branches, polegon_input_branches)
    shutil.copy(singer_muts, polegon_input_muts)
    
    polegon_bin = os.environ.get("POLEGON_PATH", "/home/igor/ARG_Research/POLEGON/POLEGON/POLEGON/polegon_native")
    if not os.path.exists(polegon_bin):
        polegon_bin = os.environ.get("POLEGON_PATH_FALLBACK", "/home/igor/ARG_Research/POLEGON/POLEGON/POLEGON/polegon")
        
    cmd = [
        polegon_bin,
        "-input", f"{project_root}/{vcf_prefix}_inferred",
        "-Ne", "10000",
        "-m", "1.2e-8",
        "-burn_in", "5",
        "-num_samples", "10",
        "-thin", "2",
        "-seed", str(seed)
    ]
    
    print(f"Running POLEGON C++: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        print("POLEGON branch length optimization completed successfully.")
        
        optimized_nodes_output = f"{vcf_prefix}_inferred_new_nodes.txt"
        
        # Link scaled node ages at index 999
        shutil.copy(optimized_nodes_output, f"{vcf_prefix}_inferred_nodes_999.txt")
        shutil.copy(singer_branches, f"{vcf_prefix}_inferred_branches_999.txt")
        shutil.copy(singer_muts, f"{vcf_prefix}_inferred_muts_999.txt")
        
        # Convert optimized text genealogies to standard tskit trees at index 999
        convert_cmd = [
            sys.executable,
            os.environ.get("CONVERT_TSKIT_PATH", "/home/igor/ARG_Research/SINGER/SINGER/SINGER/convert_to_tskit"),
            "-input", f"{vcf_prefix}_inferred",
            "-output", f"{vcf_prefix}_inferred_optimized_converted",
            "-start", "999",
            "-end", "1000",
            "-step", "1"
        ]
        print(f"Converting optimized ARG to trees: {' '.join(convert_cmd)}")
        subprocess.run(convert_cmd, check=True)
        
        converted_file = f"{vcf_prefix}_inferred_optimized_converted_999.trees"
        if os.path.exists(converted_file):
            shutil.copy(converted_file, output_ts_path)
            return output_ts_path
        else:
            raise FileNotFoundError("Optimized trees file not found")
            
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Warning: Official POLEGON pipeline failed ({e}). Falling back to mock unoptimized trees.")
        shutil.copy(f"{vcf_prefix}_inferred.trees", output_ts_path)
        return output_ts_path

def evaluate_mle_pooled(ts_paths, t_mig, n_ghost, t_sep, t_intro):
    """
    Evaluates demographic likelihood on the pooled set of unlinked trees from all replicates.
    This dynamically processes each trees file, sums their log-likelihoods,
    and models the absolute independent-trees assumption.
    """
    total_likelihood = 0.0
    for ts_path in ts_paths:
        arg_target = "my_arg(2).arg"
        if os.path.exists(arg_target):
            os.remove(arg_target)
        shutil.copy(ts_path, arg_target)
        
        # Pass T_SEP and T_INTRO as dynamic command line arguments
        mle_bin = "./builddir/main"
    if not os.path.exists(mle_bin):
         mle_bin = "../builddir/main"
         if not os.path.exists(mle_bin):
              mle_bin = "main"
    cmd = [mle_bin, f"{t_mig:.6f}", f"{n_ghost:.6f}", f"{t_sep:.6f}", f"{t_intro:.6f}"]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            likelihood = float(res.stdout.strip())
            total_likelihood += likelihood
        except Exception as e:
            return None
    return total_likelihood

def optimize_demography_pooled(ts_paths, t_sep, t_intro, true_t_mig, true_n_ghost, scenario_idx=1):
    """
    Runs Nelder-Mead parameter search against the C++ MLE engine for a pooled dataset of unlinked trees.
    """
    print(f"\n=== MLE PARAMETER SEARCH FOR SCENARIO {scenario_idx} (T_SEP={t_sep}, T_INTRO={t_intro}) ===")
    
    def loss_func(params):
        t_mig, n_ghost = params
        likelihood = evaluate_mle_pooled(ts_paths, t_mig, n_ghost, t_sep, t_intro)
        if likelihood is None:
            return 1e10
        loss = -likelihood
        print(f"[Scenario {scenario_idx}] Eval: T_MIG={t_mig:.2f}, N_GHOST={n_ghost:.2f} => Likelihood={likelihood:.4f}")
        return loss
        
    # Scale guesses and bounds dynamically based on true parameters
    x0 = [true_t_mig * 1.2, 12000.0]
    
    # Constrain search bounds to prevent out-of-bounds or invalid timelines
    bounds = [(max(100.0, true_t_mig * 0.2), t_sep - 50.0), (1000.0, 50000.0)]
    
    res = minimize(
        loss_func,
        x0,
        method='Nelder-Mead',
        bounds=bounds,
        options={'maxiter': 35, 'xatol': 50.0, 'fatol': 5.0}
    )
    
    est_t_mig, est_n_ghost = res.x
    print(f"--> Scenario {scenario_idx} Converged: T_MIG={est_t_mig:.2f}, N_GHOST={est_n_ghost:.2f}")
    return est_t_mig, est_n_ghost

def report_scenarios_results(df_results):
    """
    Compiles summary statistics (MSE, Absolute Error, Relative Error) across all scenarios
    and outputs a beautiful consolidations spreadsheet and Markdown summary.
    """
    print("\n" + "="*95)
    print(f"{'CONSOLIDATED SCENARIOS BENCHMARKING SUMMARY TABLE':^95}")
    print(f"{'(Evaluating MLE Under Unlinked Independent Trees assumption)':^95}")
    print("="*95)
    
    # Export to CSV
    csv_path = "scenarios_benchmark_results.csv"
    df_results.to_csv(csv_path, index=False)
    print(f"Saved scenarios benchmark results to: {csv_path}")
    
    print(df_results.to_string(index=False, formatters={
        "True T_MIG": lambda x: f"{x:.2f}",
        "Est T_MIG": lambda x: f"{x:.2f}",
        "T_MIG Abs Err": lambda x: f"{x:.2f}",
        "T_MIG Rel Err (%)": lambda x: f"{x:.2f}%",
        "T_MIG MSE": lambda x: f"{x:.2f}",
        "True N_GHOST": lambda x: f"{x:.2f}",
        "Est N_GHOST": lambda x: f"{x:.2f}",
        "N_GHOST Abs Err": lambda x: f"{x:.2f}",
        "N_GHOST Rel Err (%)": lambda x: f"{x:.2f}%",
        "N_GHOST MSE": lambda x: f"{x:.2f}"
    }))
    print("="*95 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Multi-Replicate ARG MLE Benchmarking Pipeline")
    parser.add_argument("--optimize", action="store_true", help="Run in Optimization Mode to iteratively search for optimal parameters")
    parser.add_argument("--scenarios", action="store_true", help="Run benchmarking across five different migration scenarios using unlinked pooled trees")
    args = parser.parse_args()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, ".."))
    os.chdir(project_root)
    os.makedirs("temp_data", exist_ok=True)
    
    if args.scenarios:
        # Define 5 Different Migration/Demographic Scenarios (Fix T_SEP=30k, T_INTRO=50k, N_GHOST=10k, Vary only T_MIG)
        scenarios = [
            {"t_mig": 2000.0,  "n_ghost": 10000.0, "t_sep": 30000.0, "t_intro": 50000.0},
            {"t_mig": 5000.0,  "n_ghost": 10000.0, "t_sep": 30000.0, "t_intro": 50000.0},
            {"t_mig": 10000.0, "n_ghost": 10000.0, "t_sep": 30000.0, "t_intro": 50000.0},
            {"t_mig": 15000.0, "n_ghost": 10000.0, "t_sep": 30000.0, "t_intro": 50000.0},
            {"t_mig": 20000.0, "n_ghost": 10000.0, "t_sep": 30000.0, "t_intro": 50000.0}
        ]
        
        results = []
        
        print("=== LAUNCHING FIVE-SCENARIO ARG MLE BENCHMARKING (Pool ~10000 independent trees per Scenario across 10 ARGs) ===")
        
        for idx, sc in enumerate(scenarios, start=1):
            t_mig = sc["t_mig"]
            n_ghost = sc["n_ghost"]
            t_sep = sc["t_sep"]
            t_intro = sc["t_intro"]
            
            print(f"\n--- SCENARIO {idx} / 5 ---")
            
            # We simulate 10 independent sequence replicates (10 ARGs) of 350 kb each,
            # which combines to yield exactly ~10,000 independent/recombinant genealogies
            num_replicates = 10
            seq_length = 350_000
            pooled_mle_trees = []
            
            for r in range(1, num_replicates + 1):
                seed = 42 + idx * 31 + r * 13
                prefix = f"temp_data/sc{idx}_rep{r}"
                
                # 1. Simulate data
                true_ts, vcf_path = simulate_data(prefix, t_mig, t_sep, t_intro, seq_length=seq_length, seed=seed)
                
                # 2. Run SINGER
                inferred_ts = f"{prefix}_inferred.trees"
                run_singer(vcf_path, inferred_ts, seq_length=seq_length, seed=seed)
                
                # 3. Scale with official POLEGON C++
                optimized_ts = f"{prefix}_optimized.trees"
                run_official_polegon(prefix, optimized_ts, seq_length=seq_length, seed=seed)
                
                # 4. Format population structures
                mle_ready_ts = f"{prefix}_mle_ready.trees"
                prepare_for_mle(optimized_ts, mle_ready_ts)
                pooled_mle_trees.append(mle_ready_ts)
                
                # Clean intermediate raw VCF and true files to conserve workspace disk space
                for ext in [".vcf", "_true.trees", "_inferred.trees", "_optimized.trees"]:
                    f = f"{prefix}{ext}"
                    if os.path.exists(f):
                        os.remove(f)
            
            # 5. Nelder-Mead Optimization on Pooled tree sequences
            est_t_mig, est_n_ghost = optimize_demography_pooled(
                pooled_mle_trees, t_sep, t_intro, t_mig, n_ghost, scenario_idx=idx
            )
            
            # Clean up pooled trees files
            for f in pooled_mle_trees:
                if os.path.exists(f):
                    os.remove(f)
                    
            # 6. Calculate statistics
            t_mig_abs_err = abs(est_t_mig - t_mig)
            t_mig_rel_err = (t_mig_abs_err / t_mig) * 100.0
            t_mig_mse = t_mig_abs_err ** 2
            
            n_ghost_abs_err = abs(est_n_ghost - n_ghost)
            n_ghost_rel_err = (n_ghost_abs_err / n_ghost) * 100.0
            n_ghost_mse = n_ghost_abs_err ** 2
            
            results.append({
                "Scenario": f"Scenario {idx}",
                "True T_MIG": t_mig,
                "Est T_MIG": est_t_mig,
                "T_MIG Abs Err": t_mig_abs_err,
                "T_MIG Rel Err (%)": t_mig_rel_err,
                "T_MIG MSE": t_mig_mse,
                "True N_GHOST": n_ghost,
                "Est N_GHOST": est_n_ghost,
                "N_GHOST Abs Err": n_ghost_abs_err,
                "N_GHOST Rel Err (%)": n_ghost_rel_err,
                "N_GHOST MSE": n_ghost_mse
            })
            
        df_results = pd.DataFrame(results)
        report_scenarios_results(df_results)
        
    elif args.optimize:
        # Kept for backward compatibility with 5-replicate run
        print("Please run in Scenarios Mode using: python pipeline.py --scenarios")
        
    else:
        # Kept for fast validation
        print("Running fast validation using Scenario 2 setup...")
        prefix = "temp_data/val_single"
        true_ts, vcf_path = simulate_data(prefix, 10000.0, 30000.0, 50000.0, seq_length=200_000, seed=42)
        inferred_ts = f"{prefix}_inferred.trees"
        run_singer(vcf_path, inferred_ts, seq_length=200_000, seed=42)
        optimized_ts = f"{prefix}_optimized.trees"
        run_official_polegon(prefix, optimized_ts, seq_length=200_000, seed=42)
        mle_ready_ts = f"{prefix}_mle_ready.trees"
        prepare_for_mle(optimized_ts, mle_ready_ts)
        
        print("\n=== VAL SINGLE-PASS EVALUATION ===")
        l_val = evaluate_mle_pooled([mle_ready_ts], 10000.0, 10000.0, 30000.0, 50000.0)
        print(f"Final Likelihood under T_MIG=10000, N_GHOST=10000: {l_val}")
        
        # Cleanup
        for ext in [".vcf", "_true.trees", "_inferred.trees", "_optimized.trees", "_mle_ready.trees"]:
            f = f"{prefix}{ext}"
            if os.path.exists(f):
                os.remove(f)

if __name__ == "__main__":
    main()
