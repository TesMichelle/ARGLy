import sys
import os
import subprocess
import tskit
import msprime
import pandas as pd
import numpy as np
import shutil
from scipy.optimize import minimize

# Target demographic parameters
TRUE_T_MIG = 10000.0
TRUE_N_GHOST = 10000.0
FIXED_T_SEP = 30000.0
FIXED_T_INTRO = 50000.0

def simulate_tree(seed):
    """
    Simulates a single tree (recombination = 0) of length 2 Mb.
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
        samples={"AFR": 25, "SAM": 25},  # 50 genomes total
        demography=demography,
        sequence_length=2_000_000,
        recombination_rate=0.0,
        random_seed=seed
    )
    return ts

def combine_tree_sequences(ts_list):
    combined_tables = tskit.TableCollection(sequence_length=0.0)
    combined_tables.populations.add_row()
    
    node_offset = 0
    total_length = 0.0
    
    for ts in ts_list:
        tables = ts.dump_tables()
        L = ts.sequence_length
        
        # Append nodes, dropping individual reference and mapping population
        for node in tables.nodes:
            orig_pop = node.population
            if orig_pop == 1:  # SAM
                pop_id = tskit.NULL
            else:
                pop_id = 0
                
            combined_tables.nodes.add_row(
                flags=node.flags,
                time=node.time,
                population=pop_id,
                individual=tskit.NULL,
                metadata=node.metadata
            )
            
        # Append edges
        for edge in tables.edges:
            combined_tables.edges.add_row(
                left=edge.left + total_length,
                right=edge.right + total_length,
                parent=edge.parent + node_offset,
                child=edge.child + node_offset,
                metadata=edge.metadata
            )
            
        node_offset += len(tables.nodes)
        total_length += L
        
    combined_tables.sequence_length = total_length
    combined_tables.sort()
    
    return combined_tables.tree_sequence()

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

def optimize_replicate(ts_path):
    def loss_func(params):
        t_mig, n_ghost = params
        likelihood = evaluate_mle(ts_path, t_mig, n_ghost)
        return -likelihood if likelihood is not None else 1e10
        
    x0 = [12000.0, 12000.0]
    bounds = [(100.0, FIXED_T_SEP - 50.0), (1000.0, 50000.0)]
    res = minimize(loss_func, x0, method='Nelder-Mead', bounds=bounds, options={'maxiter': 50, 'xatol': 50.0, 'fatol': 5.0})
    return res.x

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, ".."))
    os.chdir(project_root)
    num_replicates = 100
    num_trees_per_rep = 50
    estimates = []
    
    print(f"=== LAUNCHING OPTIMIZED TRUE TREES BOOTSTRAP ({num_replicates} Replicates) ===")
    
    for r in range(1, num_replicates + 1):
        try:
            # 1. Simulate 50 independent tree sequences
            ts_list = []
            for t_idx in range(1, num_trees_per_rep + 1):
                seed = r * 1000 + t_idx * 17 + 99
                ts = simulate_tree(seed)
                ts_list.append(ts)
                
            # 2. Combine into a single combined tree sequence file
            combined_ts = combine_tree_sequences(ts_list)
            combined_ts_path = f"boot_true_rep{r}_combined.trees"
            combined_ts.dump(combined_ts_path)
            
            # 3. Optimize parameters on the combined tree sequence
            est_t_mig, est_n_ghost = optimize_replicate(combined_ts_path)
            estimates.append({"Replicate": r, "Est_T_MIG": est_t_mig, "Est_N_GHOST": est_n_ghost})
            print(f"[Rep {r}/{num_replicates}] Completed: Est T_MIG={est_t_mig:.2f}, Est N_GHOST={est_n_ghost:.2f}")
            
            # Clean up the combined file
            if os.path.exists(combined_ts_path):
                os.remove(combined_ts_path)
                
        except Exception as e:
            print(f"[Rep {r}/{num_replicates}] Failed with error: {e}")
            
    if os.path.exists("my_arg(2).arg"):
        os.remove("my_arg(2).arg")
        
    df = pd.DataFrame(estimates)
    df.to_csv("bootstrap_true_estimates.csv", index=False)
    
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
    print(f"{'BOOTSTRAP SUMMARY STATISTICS (B = 100, 50 True Trees)':^80}")
    print("="*80)
    print(f"T_MIG:   Mean = {mean_t_mig:.2f}, SE = {se_t_mig:.2f}")
    print(f"         95% CI (Percentile) = [{ci_lower_t_mig:.2f}, {ci_upper_t_mig:.2f}]")
    print("-"*80)
    print(f"N_GHOST: Mean = {mean_n_ghost:.2f}, SE = {se_n_ghost:.2f}")
    print(f"         95% CI (Percentile) = [{ci_lower_n_ghost:.2f}, {ci_upper_n_ghost:.2f}]")
    print("="*80)

if __name__ == "__main__":
    main()
