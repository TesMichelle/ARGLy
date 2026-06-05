import tskit
import numpy as np
from scipy.optimize import minimize, LinearConstraint

def optimize_edge_lengths(ts_path, output_path, T_sep=1000.0, epsilon=1e-5):
    """
    Optimizes node times in a tree sequence to resolve illegal early coalescences
    under a population-split constraint (T_sep) without modifying the topology.
    """
    print(f"Loading tree sequence from {ts_path}...")
    ts = tskit.load(ts_path)
    
    num_nodes = ts.num_nodes
    t_orig = ts.tables.nodes.time.copy()
    
    # 1. Identify populations of samples
    # We want to identify sample nodes and their populations
    # 0=AFR (or other population), 1=SAM (or outgroup)
    samples = ts.samples()
    sample_pops = {u: ts.node(u).population for u in samples}
    
    # Track cross-population coalescence nodes
    # A node u is a cross-population coalescence if in any local tree,
    # its descendant samples include both AFR and SAM samples,
    # but none of its children have descendant samples of both populations.
    cross_coal_nodes = set()
    
    print("Identifying cross-population coalescence nodes...")
    for tree in ts.trees():
        # Precompute sample populations under each node in this tree
        # We can do a postorder traversal to find descendant sample sets
        descendants = {u: set() for u in tree.nodes()}
        
        # Initialize leaves
        for u in samples:
            if u in descendants:
                descendants[u].add(sample_pops[u])
                
        # Postorder traversal to propagate populations up
        for u in tree.postorder():
            for child in tree.children(u):
                descendants[u].update(descendants[child])
            
            # Check if this node is a cross-population coalescence
            # Must have both populations (0 and 1, or generally >= 2 distinct populations among descendants)
            if len(descendants[u]) >= 2:
                # Check if any child already has both
                child_has_both = False
                for child in tree.children(u):
                    if len(descendants[child]) >= 2:
                        child_has_both = True
                        break
                if not child_has_both:
                    cross_coal_nodes.add(u)
                    
    print(f"Found {len(cross_coal_nodes)} cross-population coalescence nodes.")
    
    # 2. Build Optimization Matrices
    # Variables: x = new node times (size = num_nodes)
    # Objective: minimize sum((x - t_orig)^2)
    # Bounds: for samples, bound is (0, 0); for others, (0, None)
    bounds = []
    for u in range(num_nodes):
        node = ts.node(u)
        if node.is_sample():
            bounds.append((0.0, 0.0))
        else:
            bounds.append((0.0, None))
            
    # Constraints: A * x >= b
    # A_rows will contain the constraint equations
    A_list = []
    b_list = []
    
    # A. Edge constraints (parent - child >= epsilon)
    # Collect unique parent-child relationships from edges table
    edges = ts.tables.edges
    unique_edges = set(zip(edges.parent, edges.child))
    
    for parent, child in unique_edges:
        # parent - child >= epsilon
        row = np.zeros(num_nodes)
        row[parent] = 1.0
        row[child] = -1.0
        A_list.append(row)
        b_list.append(epsilon)
        
    # B. Population-split constraints (u >= T_sep) for all cross_coal_nodes
    for u in cross_coal_nodes:
        row = np.zeros(num_nodes)
        row[u] = 1.0
        A_list.append(row)
        b_list.append(T_sep)
        
    # Combine into a single matrix constraint if we have constraints
    if A_list:
        A = np.vstack(A_list)
        b = np.array(b_list)
        linear_constraint = LinearConstraint(A, b, np.inf)
        constraints = [linear_constraint]
    else:
        constraints = []
        
    # 3. Solve Constrained Optimization
    print("Running SLSQP optimizer...")
    def objective(x):
        return np.sum((x - t_orig) ** 2)
        
    def jacobian(x):
        return 2.0 * (x - t_orig)
        
    # Initial guess: start from t_orig, but push cross coal nodes to at least T_sep
    x0 = t_orig.copy()
    for u in cross_coal_nodes:
        if x0[u] < T_sep:
            x0[u] = T_sep
            
    # Fix any parent-child ordering issues in initial guess
    # (Do a simple topological pass to ensure parents > children)
    tables = ts.tables
    nodes_order = ts.tables.nodes.time.argsort()
    for u in nodes_order:
        node = ts.node(u)
        if not node.is_sample():
            # parent time must be strictly greater than max child time + epsilon
            children_times = [x0[edge.child] for edge in ts.edges() if edge.parent == u]
            if children_times:
                min_parent_time = max(children_times) + epsilon
                if x0[u] < min_parent_time:
                    x0[u] = min_parent_time
                    
    res = minimize(
        objective,
        x0,
        method='SLSQP',
        jac=jacobian,
        bounds=bounds,
        constraints=constraints,
        options={'maxiter': 500, 'disp': True}
    )
    
    if not res.success:
        print(f"Optimization warning: {res.message}")
        
    # 4. Save the adjusted tree sequence
    print("Updating node times in tree sequence...")
    new_times = res.x
    
    # Create copy of tables and write new times
    new_tables = ts.dump_tables()
    new_tables.nodes.time = new_times
    
    # Sort tables to maintain tskit validity
    new_tables.sort()
    new_ts = new_tables.tree_sequence()
    
    print(f"Saving optimized tree sequence to {output_path}...")
    new_ts.dump(output_path)
    print("Done!")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python polegon.py <input.trees> <output.trees> [T_sep]")
        sys.exit(1)
    
    ts_in = sys.argv[1]
    ts_out = sys.argv[2]
    t_sep = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
    
    optimize_edge_lengths(ts_in, ts_out, T_sep=t_sep)
