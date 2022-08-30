import parse_benchmark_files as parser
import sys
import numpy as np
from time import time
import EM_PRIMA as EM
from tqdm import tqdm, trange

def main():
    tt= time()
    if len(sys.argv)>2:
        design = sys.argv[1]
        discretization = float(sys.argv[2])
    else:
        print("Please provide design name and discretization")
        exit()

    spice_file = "../benchmarks/%s.spice"%design
    spice_current_file = "../benchmarks/%s.out"%design

# tt= time()
# design = "ibmpg2"
# discretization = 10e-6
# spice_file = "%s.spice"%design
# spice_current_file = "%s.out"%design

    st1= time()
    edges, pg_unit = parser.parse_benchmark(design, spice_file, spice_current_file)
    et1 = time()


    st2 = et1
    #    print("parsed edges:",len(edges))
    seg_no, wire_no, layer_graph = parser.mapping_to_data_structure(discretization, edges, pg_unit)
    et2 = time()


    st3 = et2
    #    print("nodes:",sum([len(list(x.keys())) for x in layer_graph.values()]))
    graphs_all, count_end_nodes_all, trees_all, end_nodes = parser.create_disconnected_graphs(layer_graph)
    et3 = time()
    #print(trees_all)
    #for t in range(len(trees_all)):
    #  if trees_all[t]:
    #    if '403_551' in graphs_all[t]:
    #      print("graph")
    #      for n,node in graphs_all[t].items():
    #        print("%8s : l: %8s, r: %8s, t: %8s, b: %8s"%(n,node.get('left',''),
    #                                                      node.get('right',''),
    #                                                      node.get('top',''),
    #                                                      node.get('bottom','')))
    #        
        #if len(graphs_all[t])<500:
        #  print(graphs_all[t])

    #print("done")

    BFS_time = 0
    Matrix_Building_time = 0
    Moment_Computation_time = 0
    Krylov_time = 0
    tot_num_nodes = 0
    tot_num_wires = 0
    tot_num_seg = 0
    err_count = 0
    err = []
    relaunched = []

    for graph_no in trange(len(graphs_all)):
    #for graph_no in range(0,1):
        #print(graph_no)

        nwire = count_end_nodes_all[graph_no] - 1
        nseg = len(graphs_all[graph_no]) - 1

        graph = graphs_all[graph_no]
        tot_num_nodes  += len(graph)
        tot_num_wires  += nwire
        tot_num_seg   += nseg
        begin = time()

        start_node = None

        # search for a start_node (with no input edge from left or bottom)
        for key in graph:
            if 'left' not in graph[key] and 'bottom' not in graph[key]:
                start_node = key
                break
        #print("length:", len(graph))
        topo_sort = EM.BFS(graph, start_node)
        end = time()
        BFS_time = BFS_time + end - begin

        scale = 1e9

        begin = time()
        G, C, L, J = EM.build_matrix(nseg, nwire, graph, scale)
        end = time()

        #if len(graph) == 3:
        #  print("G mat")
        #  print(G)
        Matrix_Building_time = Matrix_Building_time + end - begin

        order = 4
        order = min (4, len(graph)-1)
        begin = time()
        Allmoments = EM.calculate_moments(nseg, graph, order, C, topo_sort)
        end = time()
        Moment_Computation_time = Moment_Computation_time + end - begin

        begin = time()
        D, r, err_flag, relaunch = EM.Krylov_subspace(nseg, nwire, order, Allmoments, G,
                                            C, L, J, len(graph), graph, topo_sort)
        if relaunch:
          relaunched.append(graph_no)

        if err_flag == 1:
            err_count = err_count + 1
            err.append(graph_no)
        end = time()
        Krylov_time = Krylov_time + end - begin

    print("Time to parse: %d"%(et1-st1))
    print("Time to create data structure: %d"%(et2-st2))
    print("Time to traverse disconnected: %d"%(et3-st3))

    print("Time for BFS: %d"%BFS_time)
    print("Time to Matrix Buildup: %d"%Matrix_Building_time)
    print("Time to Compute Moments: %d"%Moment_Computation_time)
    print("Time to Compute Pole-Zero: %d"%Krylov_time)
    print("Total Time: %d"%(end-tt))

    print("Number graphs: ",len(graphs_all))
    print("Total Number nodes: ",tot_num_nodes)
    print ("Total Number of segments", tot_num_seg)
    print ("Total Number of wires", tot_num_wires)
    print("Number errored: ",err_count)
    print("Number relaunched: ",len(relaunched))

    # printing test-graph with warning/errors
    # index = 406
    # testgraph = graphs_all[index]
    # filename = design+"_tree_"+str(index)+".txt"
    # with open(filename, 'w') as f:
    #     for key, value in testgraph.items():
    #         row = "'"+str(key)+"'"+":" + str(value) + ","
    #         f.write(row)
    #         f.write('\n')



    # # Constants
    # B = 28e9;                           # Pa = N/m^2
    # Omega = 1.18e-29;                   # m^3
    # kT = 1.38e-23*378;                  # J
    # Da = 1.3e-9*np.exp(-0.8*1.6e-19/kT);   # m^2/s
    # Z = 1;                              # dimensionless
    # q = 1.6e-19;                        # Coulomb
    # rho = 2.25e-8;                      # ohm m

    # kappa = Da*B*Omega/kT;
    # beta = Z*q*rho/Omega;

    # # Computing stress using pole-zero (D,r)

    # timearray = np.array([1e7, 1e8, 6.31e8])
    # stress = np.zeros((nwire+1, len(timearray)))

    # for n in range(nwire+1):
    #     for m in range(order):
    #         if D[m] != 0:
    #             stress[n,:] = stress[n,:] + r[n,m]*(1 - np.exp(-timearray*kappa*scale/D[m]))


main()
