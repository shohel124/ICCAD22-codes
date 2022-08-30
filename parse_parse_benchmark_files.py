import numpy as np
import warnings
warnings.filterwarnings("error")

def build_matrix(nseg, nwire, graph, scale):
# Defining G, C, L, J matrices

    G = np.zeros((nseg+1, nseg+1))
    C = np.zeros(nseg+1) # replacement of 2D diagonal matrix
    WireEndPoints = np.zeros(nwire+1)
    L = np.zeros((nseg+1,nwire+1))
    J = np.zeros(nseg+1)

#    scale = 1e9 #C matrix scale factor
    for key in graph:
        i = graph[key]['id']

        if 'right' in graph[key]:
            k = graph[graph[key]['right']]['id']
            dx = graph[key]['right_len']
            Gval = 1/dx
            Cval = 0.5*scale*dx
            G[i,i] = G[i,i] + Gval
            G[k,k] = G[k,k] + Gval
            G[i,k] = G[i,k] - Gval
            G[k,i] = G[k,i] - Gval
            C[i] = C[i] + Cval
            C[k] = C[k] + Cval

        if 'top' in graph[key]:
            k = graph[graph[key]['top']]['id']
            dx = graph[key]['top_len']
            Gval = 1/dx
            Cval = 0.5*scale*dx
            G[i,i] = G[i,i] + Gval
            G[k,k] = G[k,k] + Gval
            G[i,k] = G[i,k] - Gval
            G[k,i] = G[k,i] - Gval
            C[i] = C[i] + Cval
            C[k] = C[k] + Cval

    for key in graph:
        if graph[key]['wire_end'] == 1:
            k = int(graph[key]['id'])
            i = int(graph[key]['wire_end_id'])
            WireEndPoints[i] = k
            L[k,i] = 1
            if 'right' in graph[key]:
                J[k] = J[k] + graph[key]['right_cur']
            if 'left' in graph[key]:
                J[k] = J[k] + graph[key]['left_cur']
            if 'bottom' in graph[key]:
                J[k] = J[k] + graph[key]['bottom_cur']
            if 'top' in graph[key]:
                J[k] = J[k] + graph[key]['top_cur']

# Replacing the last row in G, C matrices to avoid singularity
##C[-1] = 0 # last row of C do not need to be forced to zero, skipped during iterations in matrix operations
    G[-1,:] = C
    J[-1]=0 #replacing last row in J with 0

    return G,C,L,J

def BFS(graph, start_node):

    # Mark all the vertices as not visited
    node_no = len(graph)
    visited = [False]*node_no
    queue = []
    BFsort = []

    queue.append(start_node) # (starting node)
    start_node_id = int(graph[start_node]['id'])

    visited[start_node_id] = True

    while queue:
        v = queue.pop(0)
        BFsort.append(v)

        if 'right' in graph[v]:
            i = graph[graph[v]['right']]['id']
            if visited[i] == False:
                queue.append(graph[v]['right'])
                visited[i] = True

        if 'left' in graph[v]:
            i = graph[graph[v]['left']]['id']
            if visited[i] == False:
                queue.append(graph[v]['left'])
                visited[i] = True

        if 'bottom' in graph[v]:
            i = graph[graph[v]['bottom']]['id']
            if visited[i] == False:
                queue.append(graph[v]['bottom'])
                visited[i] = True

        if 'top' in graph[v]:
            i = graph[graph[v]['top']]['id']
            if visited[i] == False:
                queue.append(graph[v]['top'])
                visited[i] = True

    return BFsort


def calculate_moments(nseg, graph, order, C, topo_sort):

    moment = np.zeros((order, nseg+1))
    for l in range(order) :
        j = np.zeros(nseg)
        offset = np.zeros(nseg+1)
        jj = np.zeros(nseg+1)

        if l!=0:  #for computing m1, m2, m3
            j = np.multiply(C,moment[l-1,:])
    #        print(j)

            jj = np.zeros(nseg+1)
            for i in reversed(range(1,nseg + 1)):
                k = topo_sort[i]
                index = int(graph[k]['id'])
    #            print(index, k)

                if 'right' not in graph[k] and 'bottom' not in graph[k]:
                    jj[index] = j[index]

                if 'right' in graph[k] and 'bottom' not in graph[k]:
                    right_index = int(graph[graph[k]['right']]['id'])
                    jj[index] = j[index] + jj[right_index]

                if 'right' not in graph[k] and 'bottom' in graph[k]:
                    bottom_index = int(graph[graph[k]['bottom']]['id'])
                    jj[index] = j[index] + jj[bottom_index]

                if 'right' in graph[k] and 'bottom' in graph[k]:
                    right_index = int(graph[graph[k]['right']]['id'])
                    bottom_index = int(graph[graph[k]['bottom']]['id'])
                    jj[index] = j[index] + jj[right_index] + jj[bottom_index]

        for i in range(nseg):
            k = topo_sort[i+1]
            index = int(graph[k]['id'])


            if 'left' in graph[k]:

                left_index = int(graph[graph[k]['left']]['id'])
                length = graph[graph[k]['left']]['right_len'] #L

                if l==0:
                    cur_den = graph[graph[k]['left']]['right_cur'] #beta*j
                    offset[index] = offset[left_index] - cur_den*length

                else:
                    offset[index] = offset[left_index] - jj[index]*length

            if 'top' in graph[k] :

                top_index = int(graph[graph[k]['top']]['id'])
                length = graph[graph[k]['top']]['bottom_len'] #L

                if l==0:
                    cur_den = graph[graph[k]['top']]['bottom_cur'] #beta*j
                    offset[index] = offset[top_index] - cur_den*length

                else:
                    offset[index] = offset[top_index] - jj[index]*length

        RHS = -np.dot(C, offset)
        LHS = np.sum(C)
        moment[l,0]=RHS/LHS

        moment[l,:] = moment[l,0] + offset
       # if len(graph) ==3:
       #   print("Order:",l)
       #   print("J:\n",j)
       #   print("jj:\n",jj)
       #   print("moment:\n",moment)
       #   print("offset:\n",offset)
       #   print("#######")
    Allmoments = moment.transpose()

    return Allmoments

def Krylov_subspace(nseg, nwire, order, Allmoments, G, C, L, J, g_len,graph,topo_sort):
    err_flag = 0
    Xall = np.zeros((nseg+1,order))
    Xall[:,0] = Allmoments[:,0]
    D = np.zeros(order)
    r = np.zeros((nwire+1, order))
    relaunch = 0
    try:
      for l in range(order):
          if l > 0:
              v=Allmoments[:,l]
              for k in range(l):
                  num = np.dot(v,Allmoments[:,k])
                  denom = np.dot(Allmoments[:,k],Allmoments[:,k])
                  v = v - num/denom * Allmoments[:,k]
              Xall[:,l] = v

      X = Xall
      Xt = np.transpose(X)
      P = np.matmul(G,X)
      Gnew = np.matmul(Xt,P)
      #if g_len == 3:
      #  print("modified G mat")
      #  print(Gnew)
      #  print("X mat")
      #  print(Xall)
      #Q = np.matmul(C,X)
      #Since C is diagonal, Q can be calculated as Q[ij] = C[ii] X[ij]
      Q = np.zeros((nseg+1,order))
      Cnew = np.zeros((order,order))

      # for i in range(nseg): # iterate 1 time less as the final row(element) of C matrix(array) is zero
      #     Q[i,:]=C[i]*X[i,:]

      for i in range(order):
          Q[0:nseg,i]=np.multiply(C[0:nseg],Xt[i,0:nseg])

      for i in range(order):
          for k in range(order):
              Cnew[i,k] = np.dot(Xt[i,:],Q[:,k])

      #Cnew = np.matmul(np.transpose(X),Q)
      #Since C is diagonal, Cnew can be calculated as Cnew[ij] = C[kk] X[ki] X[kj]

      # for i in range(order):
      #     for j in range(order):
      #         for k in range(nseg): # iterate 1 time less as the final row(element) of C matrix(array) is zero
      #             Cnew[i,j] = Cnew[i,j] + C[k]*X[k,i]*X[k,j]

      #print("computing inverse G")
      inv_G = np.linalg.inv(Gnew)
      #print("InvG completed")
      A = np.matmul(inv_G, Cnew)
      #A = np.matmul(np.linalg.inv(Gnew),Cnew)

      D, S = np.linalg.eig(A)

      # eigen-value decomposition yields A = S * D * inv(S)
      # D is a diagonal matrix in linear algebra
      # For efficient storage, eig function in Python returns a 1-D array with the diagonal elements
      # The diagonal elements also represents the poles

      Lnew = np.matmul(Xt,L)
      Bnew = np.matmul(Xt,J)
      #print("computing inverse S")
      inv_S = np.linalg.inv(S)
      #print("InvS completed")
      nu = np.matmul (inv_S,inv_G)
      #nu = np.matmul(np.linalg.inv(S), np.linalg.inv(Gnew))
      my_nu = np.matmul(nu,Bnew)

      for k in range(nwire+1):
          my_mu = np.matmul(np.transpose(S),Lnew[:,k])
          for l in range(order):
              r[k,l] = my_mu[l]*np.transpose(my_nu[l])

    except Exception as e:
        #print("Failing for G nodes:",g_len)
        #print("Moments:", Allmoments)
        #print("G",G)
        #print("X",X)
        #print("G_hat",G)
        relaunch = 1
        if order >1:
          order = order-1
          Allmoments = calculate_moments(nseg, graph, order, C, topo_sort)
          D, r, err_flag,_ = Krylov_subspace(nseg, nwire, order, Allmoments, G,
                                           C, L, J, len(graph), graph, topo_sort)
        else:
          print("Graph fails with order 1")
          err_flag = 1

    # return err_flag
    return D, r, err_flag, relaunch
