import re
import numpy as np
_prefix = {
 'f': 1e-15,  # femto
 'p': 1e-12,  # pico
 'n': 1e-9,   # nano
 'u': 1e-6,   # micro
 'm': 1e-3,   # milli
 }

# Widths not specified in spice file calculated manually
# width in um
ibm_benchmark_widths={
'ibmpg1': {
  'M5': 10.8, # (33+183)/20 see IBM paper
  'M6':  14 #(92+188)/20 
  },

'ibmpg2':{
  'M3': 3.6, #(72/20)
  'M4': 2.4, # 48/20
  'M5': 3.6, # 72/20
  'M6': 2.4  # 48/20 ## 7?
  },
'ibmpg3': {
  'M2': 4.32, # 86.4/20
  'M3': 6.48, #129.6/20
  'M4': 4.32, # 86.4/20
  'M5': 6.48, #129.6/20
  'M6': 4.32  # 86.4/20
  },
'ibmpg4': {
  'M1': 3.60, #72/20
  'M2': 2.40, #48/20
  'M3': 3.60, #72/20
  'M4': 2.40, #48/20
  'M5': 3.60, #72/20
  'M6': 2.40  #48/20
  },
'ibmpg5': {
  'M4': 4.10, #82
  'M5': 4.10, #82
  'M6': 4.10  #82
  },
'ibmpg6': {
  'M4': 2.01, #40.32
  'M5': 2.01, #40.32
  'M6': 2.01  #40.32
  },
'test': {
  'M5': 3.6, #40.32
  'M6': 2.4  #40.32
  },
'test_new': {
  'M3': 3.6, #(72/20)
  },
 'unit_test':{
  'M4': 1e6
 }
}

#from A High Performance 180 nm Generation Logic Technology paper
# in um

ibm_benchmark_thickness={ # Back calculated with R,w and t with rho of 2.67e-8ohm/m values
'ibmpg1': {
  'M5': 0.866, # 1.13,
  'M6': 3 
  },

'ibmpg2':{
  'M3': 0.68, #0.55,
  'M4': 0.68, 
  'M5': 1.85, 
  'M6': 17.52 
  },
'ibmpg3': {
  'M2': 0.25, 
  'M3': 0.57, 
  'M4': 1.23, 
  'M5': 1.68, 
  'M6': 30.32  # 3.032
  },
'ibmpg4': {
  'M1': 0.22, 
  'M2': 0.34, 
  'M3': 0.82, 
  'M4': 0.68, 
  'M5': 2.60, 
  'M6': 17.5  
  },
'ibmpg5': {
  'M4': 0.30, 
  'M5': 0.65, 
  'M6': 10.25  
  },
'ibmpg6': {
  'M4': 0.28, #40.32
  'M5': 0.46, #40.32
  'M6': 30  #40.32
  },
'test': {
  'M5': 3.6, #40.32
  'M6': 2.4  #40.32
  },
'test_new':{
  'M3': 0.68, #0.55,
  } ,
 'unit_test':{
  'M4': 1e6
 }
}

# Constants
B = 28e9;                           # Pa = N/m^2
Omega = 1.18e-29;                   # m^3
kT = 1.38e-23*378;                  # J
Da = 1.3e-9*np.exp(-0.8*1.6e-19/kT);   # m^2/s
Z = 1;                              # dimensionless
q = 1.6e-19;                        # Coulomb
rho = 2.25e-8;                      # ohm m

kappa = Da*B*Omega/kT;
beta = Z*q*rho/Omega;

#ibm_benchmark_thickness={  
#  'M1': 0.48,
#  'M2': 0.70,
#  'M3': 0.70,
#  'M4': 1.08,
#  'M5': 1.60,
#  'M6': 1.72
#}

def unit_convert(in_value): 
  #print("input", in_value)
  try:
    #print(in_value)
    unit = re.search('([fpnum])', in_value).group(1)
    scale = _prefix[unit]
    value  = float(in_value.replace(unit,''))*scale
    #print(value)
  except AttributeError:
    value  = float(in_value)
  #print("output", value)
  return value

def parse_benchmark(design, spice_file, spice_out_file):
  layers={}
  edges = {}
  nodes = {}
  layer_name = 'M6'
  edges = {}
  nodes= {}
  node_id = 0
  net_nodes = {}
  via_nets = []

  with open(spice_file) as fp:
    for line in fp:
      m = re.search("\s*\* layer: (.*?),(.*?) ",line)
      if m is not None:
        layer_name = m.group(1)
      m = re.search("^[rRV]",line)
      if m is not None :
        cir_type = m.group(0)
        words = line.split()
        name = words[0].lower()
        loc_str1 = words[1]
        loc_str2 = words[2]
        #TODO hack to ignore _X
        m = re.search("^(n\d+)_(\d+)_(\d+)",loc_str1)
        if m is not None:
            net_name1 =  m.group(1)
            x_loc1 = int(m.group(2))
            y_loc1 = int(m.group(3))
        else:
          via_nets.append(name)
          #print("format not found %s"%loc_str1)
          continue
        
        m = re.search("^(n\d+)_(\d+)_(\d+)",loc_str2)
        if m is not None:
            net_name2 =  m.group(1)
            x_loc2 = int(m.group(2))
            y_loc2 = int(m.group(3))
        else:
          via_nets.append(name)
          #print("format not found %s"%loc_str2)
          continue

        if (loc_str1 not in nodes) and (loc_str2 not in nodes):#both are not present
            nodes[loc_str1] = node_id
            node_id = node_id + 1
            nodes[loc_str2] = node_id
            node_id = node_id + 1
        elif loc_str1 not in nodes: #only loc1 not present
            nodes[loc_str1] = node_id
            node_id = node_id + 1
        elif loc_str2 not in nodes: #only loc2 not present
            nodes[loc_str2] = node_id
            node_id = node_id + 1
        if net_name1 not in net_nodes:
            net_nodes[net_name1] = {}
            net_nodes[net_name1]['nodes'] = {}
            net_nodes[net_name1]['connections'] = {}
        net_nodes[net_name1]['nodes'][loc_str1] = nodes[loc_str1]

        if net_name2 not in net_nodes:
            net_nodes[net_name2] = {}
            net_nodes[net_name2]['nodes'] = {}
            net_nodes[net_name2]['connections'] = {}
        net_nodes[net_name2]['nodes'][loc_str2] = nodes[loc_str2]
        dx = abs(x_loc1 - x_loc2)
        dy = abs(y_loc1 - y_loc2)
        length = dx+dy
        if  net_name1 != net_name2 :
          length = 0
          via_nets.append(name)
          continue ## do not process the vias just break and continue the loop
          net_nodes[net_name2]['connections'][net_name1] = 1
          net_nodes[net_name1]['connections'][net_name2] = 1
        elif length ==0: #via on the same net
          via_nets.append(name)
          continue
        else:
          pass
        #width = ibm_benchmark_widths[design][layer_name]
        ###################################
        # Al vs Cu
        ###################################
        #if design != 'unit_test':
        #    width = width*2.25/2.67
        ###################################
        #thickness = ibm_benchmark_thickness[design][layer_name] 
        resistance = float(words[3])

        if design == 'ibmpg3' or design == 'ibmpg6':
          pg_unit = 1e9
          length =length/1000 # convert to um from nm
        else:
          pg_unit = 1e6
          
        WT = 2.25e-8*length*1e-6/resistance
        
        if cir_type != 'V': # 0 resistance via so dont create node or edge
          edges[name] = {}
          edges[name]['net'] = net_name1
          edges[name]['node1'] = loc_str1
          edges[name]['node2'] = loc_str2
          edges[name]['loc'] = (x_loc1,y_loc1) 
          edges[name]['loc2'] = (x_loc2,y_loc2) 
          edges[name]['length'] = length * 1e-6 # um
          edges[name]['WT'] = WT
          #edges[name]['width'] = width*1e-6 # um
          #edges[name]['thickness'] = thickness*1e-6 # um
          edges[name]['res'] = resistance
          edges[name]['layer'] = layer_name
          edges[name]['current'] = 0
          edges[name]['voltage_drop'] = 0
          edges[name]['power']   = 0
  
          if loc_str1 not in node_edge_dict:
              node_edge_dict [loc_str1] = []
          node_edge_dict [loc_str1].append(name)
    
          if loc_str2 not in node_edge_dict:
              node_edge_dict [loc_str2] = []
          node_edge_dict [loc_str2].append(name)
  
  print("Completed parsing spice file, starting outputfile")
  with open(spice_out_file) as fp: 
    resistors = 0
    for line in fp: 
      if resistors == 0: 
        if re.search("\*\s+resistors",line):
          resistors = 1
        continue
      try:
        net_name = re.search('0:([rR]\w*\d*)', line).group(1)
      except AttributeError:
        continue
      net_name = net_name.lower()
      edge = edges.get(net_name)
      if edge is None :
        continue
      line = line.strip()
      data = line.split()
      current = data[1]
      voltage = data[2]
      power   = data[3]
      # Divide the current to get IR drop in 10%VDD range
      if (design == 'test' or design == 'test_new'):
        current = unit_convert(current)
      elif(design == 'ibmpg1'):
        current = unit_convert(current)/4.4 
      elif (design == 'ibmpg2'):
        current = unit_convert(current)/2.1
      else:
        current = unit_convert(current)
      voltage = unit_convert(voltage)
      power   = unit_convert(power)
      #Store current density
      #edge['current'] = -1*current/(edge['width']*edge['thickness'])
      edge['voltage_drop'] = voltage
      edge['power']   = power
  print("Completed output file")
  print("Number of via edges", len(via_nets))
  return edges, pg_unit 

def mapping_to_data_structure(discretization, edges, pg_unit):
  #nodes_end = set()
  layer_graph ={}
  #graph ={}
  wire_no = 0
  seg_no = 0
  for edge_name, edge in edges.items():
    layer = "%s_%s"%(edge['layer'],edge['net'])
    if layer  not in layer_graph:
      layer_graph[layer] ={}
    x1_edge, y1_edge = edge['loc']
    x2_edge, y2_edge = edge['loc2']
    if x1_edge > x2_edge or y1_edge > y2_edge :
      tx,ty = x2_edge,y2_edge
      x2_edge,y2_edge = x1_edge,y1_edge
      x1_edge,y1_edge = tx,ty 
      current = -edge['current']
    else:
      current = edge['current']
    x1_edge = int(x1_edge)
    x2_edge = int(x2_edge)
    y1_edge = int(y1_edge)
    y2_edge = int(y2_edge)
    x_start, y_start = x1_edge, y1_edge
    pos = 0.0
    l = int(discretization*pg_unit)
    n2=0
    while pos < int(edge['length']*pg_unit):
      if abs(x1_edge-x2_edge)>0 and abs(y1_edge-y2_edge)>0:
        print("both cant change")
        print(x1_edge,x2_edge,y1_edge,y2_edge)
        exit()
      elif abs(x1_edge-x2_edge)>0:
        x_end = x_start + l
        y_end = y_start
        dir_type = "H"
      elif abs(y1_edge-y2_edge)>0:
        y_end = y_start + l
        x_end = x_start
        dir_type = "V"
      else:
        print("Cannot be equal")
        print(x1_edge,x2_edge,y1_edge,y2_edge)
        exit()
      node1 ="%d_%d"%(x_start,y_start) 
            
      pos += l
      if pos >=int(edge['length']*pg_unit):
        node2="%d_%d"%(x2_edge,y2_edge)
        length = edge['length'] - ((pos-l)/pg_unit)
      else:
        node2 ="%d_%d"%(x_end,y_end)
        length = discretization
      x_start, y_start = x_end, y_end
      seg_no = seg_no + 1
      populate_info(layer_graph[layer], node1, node2, current, 
                    length, dir_type, wire_no)
    wire_no += 1

    ne1 = "%d_%d"%(x1_edge,y1_edge)
    ne2 = "%d_%d"%(x2_edge,y2_edge)
    if ne2 != node2:
      print("mismatch: %s %s"%(node2,ne2))
    layer_graph[layer][ne1]['wire_end'] = 1
    layer_graph[layer][ne2]['wire_end'] = 1

  
  return seg_no, wire_no, layer_graph

def populate_info(graph, node1, node2, current, length, dir_type, wire_no):
  if node1 not in graph:
      graph[node1] = {}
      graph[node1]['wire_end'] = 0
  if node2 not in graph:
      graph[node2] = {}
      graph[node2]['wire_end'] = 0
  
  if dir_type == "H":
    n1 = 'left'
    n2 = 'right'
  else:
    n2 = 'top'
    n1 = 'bottom'
  
  graph[node2]['%s'%n1]=node1
  graph[node2]['%s_len'%n1]=length 
  graph[node2]['%s_cur'%n1]=-beta*current #reversing current density direction
  graph[node2]['%s_wire'%n1]= wire_no
  
  graph[node1]['%s'%n2]=node2
  graph[node1]['%s_len'%n2]=length 
  graph[node1]['%s_cur'%n2]=beta*current
  graph[node1]['%s_wire'%n2]= wire_no
  
def create_disconnected_graphs(layer_graph):
  graphs_all = []
  n_end_nodes = []
  end_nodes = []
  trees_all = []
  for layer, graph  in layer_graph.items():
    nodes_remaining = set(graph.keys())
    print("layer: %s: %d"%(layer, len(nodes_remaining)))
    while len(nodes_remaining) >0:
      graph_branch = {}
      node = next(iter(nodes_remaining))
      nodes_branch = [node]
      nodes_queue = [node]
      nodes_remaining.remove(node)
      node_id_ = 0 
      end_id = 0
      tree = 0
      end_branch= []
      while len(nodes_queue)>0:
        node = nodes_queue.pop(0)
        graph[node]['id'] = node_id_
        node_id_ += 1
        if graph[node]['wire_end'] == 1:
          graph[node]['wire_end_id'] = end_id
          end_id += 1
          end_branch.append(node)
        if(    ('left' in graph[node] or 'right'  in graph[node])
           and ('top'  in graph[node] or 'bottom' in graph[node])):
           tree = 1
        for dir_no, direction in enumerate(['left', 'right', 'top', 'bottom']):
          if direction in graph[node]:
            opp = ['right', 'left', 'bottom', 'top']
            opp_dir = opp[dir_no]
            branch_node =  graph[node][direction]
            if branch_node in nodes_remaining:
              nodes_branch.append(branch_node) 
              nodes_queue.append(branch_node)
              nodes_remaining.remove(branch_node)
              # add a key to indicate this edge was used for traversal
              graph[node]["%s_trav"%direction] = 1
              graph[branch_node]["%s_trav"%opp_dir] = 1
            else: # delete the edges if it was not used for traversing
              if ("%s_trav"%direction not in graph[node]):
                del graph[node]["%s"%direction]
                del graph[node]["%s_cur"%direction]
                del graph[node]["%s_len"%direction]
                del graph[branch_node]["%s"%opp_dir]
                del graph[branch_node]["%s_cur"%opp_dir]
                del graph[branch_node]["%s_len"%opp_dir]
              
              
      graph_branch = {k: graph[k] for k in nodes_branch}
      print(len(nodes_remaining), len(nodes_branch), "                 " , end='\r')
      graphs_all.append(graph_branch)
      n_end_nodes.append(end_id)
      end_nodes.append(end_branch)
      trees_all.append(tree)
  return graphs_all, n_end_nodes, trees_all, end_branch

def parse_voltages(voltage_file, edges, node_edge_dict):
  print("Starting voltage file")
  with open(voltage_file) as fp:
    for line in fp:
      words = line.split()
      loc_str = words[0]
      voltage = float(words[1])
      m = re.search("^(n\d+)_(\d+)_(\d+)",loc_str)
      if m is not None:

        if loc_str in node_edge_dict:
          edge_key = node_edge_dict[loc_str]
          
          for i in range(len(edge_key)):
              edge = edge_key[i]
              if edges[edge]['node1'] == loc_str :
                  edges[edge]['v1'] = voltage
              if edges[edge]['node2'] == loc_str :
                  edges[edge]['v2'] = voltage
                  
      else: #for _x_ we ignore
        continue
  #Store current density  
  for key in edges:
      edges[key]['current'] = (edges[key]['v2'] - edges[key]['v1'])/(rho*edges[key]['length'])

  return edges
