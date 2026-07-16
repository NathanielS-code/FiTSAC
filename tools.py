# this file contains all the tools the code uses

from objects import *
from printers import *

# force free assumption
def setup_SUSY_algebra(field_num, field_list, transformation_list):
  # field_num is the field in field_list that corresponds
  # convention is that we mark free spinorial indices first, then spacetime (for ID number)
  # want to make all free indices down, reset this if transformation_list does not have this convention

  # apply DaDb
  LHS_ab = transformation_list[field_num].copy() # transformation will be an expression

  expand_all_sigmas(LHS_ab)

  for term in LHS_ab.term_list:
    for relation in term.relation_list:
      if relation.relation_type == "free1":
        temp = 0 if relation.index_type == "spacetime" else 1
        tensor = relation.tensors[0]
        index = relation.indices[0]
        tensor.index_structure[temp][index] = "down" # fixing convention
        relation.relation_type = "free2" # this is b
      elif relation.relation_type[:4] == "free": # we increase the number of each free, to make room for b
        relation.relation_type = "free" + str(int(relation.relation_type[4:])+1)
        temp = 0 if relation.index_type == "spacetime" else 1
        tensor = relation.tensors[0]
        index = relation.indices[0]
        tensor.index_structure[temp][index] = "down" # fixing convention

  # now we do the Da
  length = len(LHS_ab.term_list)
  for temp in range(length): # DON'T want length to update, don't want to hit new terms
    term = LHS_ab.term_list[temp]
    # adding a temporary tensor to be removed
    Da = Tensor([[],["down"]], None, "unassigned")
    new_relation = Free_Index(Da, 0, "spinorial", 1)
    term.tensor_list.append(Da)
    term.relation_list.append(new_relation)

    # identify the field in this term
    for tensor in term.tensor_list:
      if tensor.special_type[:5] == "field":
        field_num2 = int(tensor.special_type[6:]) - 1 # the 5th spot is saved for "B" or "F"
        field = tensor
        break # found it

    # now we build up for TERMLEVELapply_identity
    transformation_copy = transformation_list[field_num2].copy()

    old_tensors = {field, Da}
    old_relations = set() # code can handle if it is a self contracted field
    new_tensors_list = []
    new_relations_list = []
    tensor_index_dict = {}
    coeff_factor_list = []

    for i in range(len(transformation_copy.term_list)):
      term2 = transformation_copy.term_list[i]
      new_tensors_list.append(term2.tensor_list)

      space_list = []
      spin_list = []

      for relation in term2.relation_list: # just storing which free's are space and spin
        if relation.relation_type[:4] == "free":
          if relation.index_type == "spacetime":
            space_list.append(relation.relation_type)
          if relation.index_type == "spinorial" and relation.relation_type != "free1":
            spin_list.append(relation.relation_type)

      space_list = sorted(space_list, key=lambda x: int(x[4:]))
      spin_list = sorted(spin_list, key=lambda x: int(x[4:]))

      new_relations = []

      # all relations are added except for free indices, these are handled with tensor_index_dict
      for relation in term2.relation_list:
        if relation.relation_type[:4] != "free":
          new_relations.append(relation)
        else: # then it is free
          tensor = relation.tensors[0]
          index = relation.indices[0]
          index_type = relation.index_type

          # we assign Da to free1 and the rest of the field to the frees
          if relation.relation_type == "free1":
            tensor_index_dict[(Da, 0, i, "spinorial")] = [tensor, index] # this one is for the Da
            tensor.index_structure[1][index] = Da.index_structure[1][0]
          else:
            lst = space_list if index_type == "spacetime" else spin_list

            from_index = lst.index(relation.relation_type)
            tensor_index_dict[(field, from_index, i, index_type)] = [tensor, index] # tuple is (tensor, index, term_num, relation.relation_type)

            j = 0 if index_type == "spacetime" else 1
            tensor.index_structure[j][index] = field.index_structure[j][from_index]

      new_relations_list.append(new_relations)

      coeff_factor_list.append(term2.coefficient)

    output_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
    LHS_ab.term_list = LHS_ab.term_list + output_terms

  trusty_reduce_eta_and_C(LHS_ab)

  # using fierz if necessary
  if field_list[field_num].special_type[5] == "F":

    trusty_use_Fierz(LHS_ab)

    LHS = LHS_ab # the coefficeints of use_Fierz make it so we don't have to double this


  else:
    # apply DbDa
    LHS_ba = LHS_ab.copy() # it is the same as LHS_ab with a b switched
    for term in LHS_ba.term_list:
      for relation in term.relation_list:
        if relation.relation_type == "free1":
          relation.relation_type = "free2"
        elif relation.relation_type == "free2":
          relation.relation_type = "free1"

    # combine into LHS
    LHS = Expression(LHS_ab.term_list + LHS_ba.term_list)

  # set up RHS
  gamma = Gamma_mu([["up"],["down","down"]])
  partial = Partial([["down"],[]])
  field = field_list[field_num].copy()
  relation = Contraction([gamma, partial], [0,0], "spacetime")
  free1 = Free_Index(gamma, 0, "spinorial", 1)
  free2 = Free_Index(gamma, 1, "spinorial", 2)

  tensor_list = [gamma, partial, field]
  relation_list = [relation,free1,free2]

  free_num = 3 # will index, adding rest of free indices now
  for i in range(len(field.index_structure[1])): # spinorial first
    new_free = Free_Index(field, i, "spinorial", free_num)
    relation_list.append(new_free)
    free_num += 1
  for i in range(len(field.index_structure[0])): # spacetime second
    new_free = Free_Index(field, i, "spacetime", free_num)
    relation_list.append(new_free)
    free_num += 1

  term = Term(2*1j, tensor_list, relation_list) # 2 is the coefficient of RHS

  RHS = Expression([term])

  return LHS, RHS


def combine_like_terms(expression, print_level=0, mode="heavy"):

  tag_list = tagmaker(expression, print_level, mode) # print_level tells it to print more statements
  # actually a list of list of the tags

  j=0
  k=0

  while j < len(expression.term_list):
    term1 = expression.term_list[j]
    k = j + 1
    while k < len(expression.term_list): # checking pairwise, don't have to check ones that have already been checked though
      term2 = expression.term_list[k]
      if tag_list[j] == tag_list[k]:

        summed_coefficient = term1.coefficient + term2.coefficient

        term1.coefficient = summed_coefficient

        del tag_list[k] # updating tag_list to match term_list
        del expression.term_list[k]

        k -= 1 # check the one at this spot again (will be new)

        if summed_coefficient == 0:
          del tag_list[j] # know j<k
          del expression.term_list[j]
          j -= 1
          break # want to go back to check this again, since there will be a new one here
      k += 1
    j += 1

def tagmaker(expression, print_level=0, mode="heavy"):
  def build_path_level(lst, graph, order, field, max_depth, current_depth=0, path=None): # relies on free field assumption

    if path is None:
      path = []

    if current_depth < max_depth:
      lstA = lst[-1]
      all_paths_built = True # default

      for objA in lstA:
        if isinstance(objA, list): # it'll be a string if ended path, won't continue
          path_so_far = list(path)
          path_so_far.append(lst[0]) # remembering the path of tensors taken to get here
          these_paths_built = build_path_level(objA, graph, order, field, max_depth, current_depth+1, path_so_far)
          if these_paths_built == False:
            all_paths_built = False
      return all_paths_built # this is to get the 1st call of build_path_level to be True or False, to evaluate whether we should keep running
    else:
      tensor = lst[0]
      lst_to_edit = lst[-1]
      all_paths_built = True # default

      path.append(tensor)

      # if this tensor has free index we immediately end
      if order.get(tensor) is not None: # getting whether it is a free or not
        if order[tensor][:6] == r"['free": # this is how the start of something with a free will look in order
          free_tag = order[tensor]

          lst_to_edit.append(free_tag) # so that it won't continue down this path

          return all_paths_built # won't check for traceback

      # not free, so we build in the direction of each relation (or give traceback)

      for edge in graph[tensor]:
        # traceback: (only triggers if not free)
        traced = False

        for i in range(len(path)):
          if edge[0] is path[i]:
            traceback_num = len(path)-1-i
            index_type = edge[1]
            lst_to_edit.append(index_type[2] + str(traceback_num))

            traced = True

            # edge[0] can only appear in path[i] once, break is justified
            break # won't run code below

        if traced == True: # then won't run code below
          continue

        # build: (only triggers if not free and not traceback)
        all_paths_built = False # this is for return statement
        elem3 = "c" # "c" for contraction
        lst_to_edit.append([edge[0], edge[1], elem3, []]) # edge[0] is the new lst[0], this whole list is objA

      return all_paths_built

  # --------------------------------------------------------------------------------------

  # this is one is for converting tensors and deep_sorting

  def convert_and_deepsort(lst):
    if isinstance(lst[0], Tensor):
      lst[0] = lst[0].special_type # convert
    lstA = lst[-1]
    for i in range(len(lstA)):
      if isinstance(lstA[i], list): # not going down ended paths
        lstA[i] = convert_and_deepsort(lstA[i]) # deepsort
    lst[-1] = sorted(lstA, key=lambda x: (isinstance(x, list), str(x)))
    return sorted(lst, key=lambda x: (isinstance(x, list), str(x)))

  # --------------------------------------------------------------------------------------

  # used for INTERGRAPH

  # find all non-self-intersecting paths that do not contain the same tensor, "largest local disjoint graph". tensor1 gets priority on graph building
  # output all the tensors in tensor2's paths in a list
  def find_local_graph(graph, tensor1, tensor2):

    end = False

    local_graph_tensors1 = {tensor1}
    local_graph_tensors2 = {tensor2}

    while end == False:
      end = True
      # now we expand. only add tensors to local_graph_tensors if they are not already in there
      for tensor in set(local_graph_tensors1): # only need to check last added tensors, but intergraph should be rare that a bit of slowness is fine
        for edge in graph[tensor]:
          if edge[0] not in local_graph_tensors1:
            local_graph_tensors1.add(edge[0])
      for tensor in set(local_graph_tensors2): # creating another set so that we don't modify original during for loop
        for edge in graph[tensor]:
          if edge[0] in local_graph_tensors1:
            local_graph_tensors1.remove(edge[0]) # remove it from there so that this remains a dead zone
          elif edge[0] not in local_graph_tensors2:
            local_graph_tensors2.add(edge[0])
            end = False # if anything got added, we continue checking

    return local_graph_tensors2

  # --------------------------------------------------------------------------------------

  if print_level == 0.5 or print_level >= 2:
    print("before tagmaker:")
    print("0 = " + expression_to_LaTeX(expression))

  tag_list = []

  terms_hit = 0
  while terms_hit < len(expression.term_list): # len will reevaluate every time

    term = expression.term_list[terms_hit]
    terms_hit += 1

    pre_length = len(expression.term_list)

    order = {}
    graph_list = []

    # 1st algorithm:
    if mode != "super heavy":
      # term pre-processing
      if print_level == 0.5 or print_level >= 2:
        print("term initial")
        print(term_to_LaTeX(term))

      TERMLEVELtrusty_reduce_eta_and_C(term)
      if TERMLEVELtrusty_auto_kill_terms(term): # if this is true we autokill
        expression.term_list.remove(term)
        terms_hit -= 1 # deindexing because we removed one
        continue # go to next term
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C and autokill")
        print(term_to_LaTeX(term))

      if mode == "heavy":
        trusty_Levi_Civita_gamma_converter(expression, term)

        if print_level == 0.5 or print_level >= 2:
          print("Levi Civita gamma converter")
          print(term_to_LaTeX(term))

      trusty_reduce_gamma_squared(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("reduced gamma^2")
        print(term_to_LaTeX(term))

      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))

      output_terms, zeroed = TERMLEVELtrusty_gamma_trace(term)
      expression.term_list += output_terms
      if zeroed: # if this is true we autokill
        expression.term_list.remove(term)
        terms_hit -= 1 # deindexing because we removed one
        continue # go to next term
      if print_level == 0.5 or print_level >= 2:
        print("reduce gamma trace")
        print(term_to_LaTeX(term))

      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))

      reduce_two_gamma_two_partial(expression, term)

      if print_level == 0.5 or print_level >= 2:
        print("reduced two gamma two partial")
        print(term_to_LaTeX(term))

      TERMLEVELtrusty_reduce_eta_and_C(term)


      if TERMLEVELtrusty_auto_kill_terms(term): # if this is true we autokill
        expression.term_list.remove(term)
        terms_hit -= 1 # deindexing because we removed one
        continue # go to next term
    else: # 2nd algorithm:

      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))
      
      trusty_reduce_gamma_squared(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("reduced gamma^2")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))

      TERMLEVELexpand_all_gamma5(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("expand gamma5")
        print(term_to_LaTeX(term))
      
      reduce_Levi_Civita_products(expression)
      if print_level == 0.5 or print_level >= 2:
        print("reduce Levi Civita products")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))
      
      trusty_reduce_long_chains(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce long chains")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))
      
      if TERMLEVELtrusty_auto_kill_terms(term): # if this is true we autokill
        expression.term_list.remove(term)
        terms_hit -= 1 # deindexing because we removed one
        continue # go to next term
      if print_level == 0.5 or print_level >= 2:
        print("autokilled term if needed")
        print(term_to_LaTeX(term))
      
      trusty_reduce_gamma_squared(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("reduced gamma^2")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))
      
      output_terms, zeroed = TERMLEVELtrusty_gamma_trace(term)
      expression.term_list += output_terms
      if zeroed: # if this is true we autokill
        expression.term_list.remove(term)
        terms_hit -= 1 # deindexing because we removed one
        continue # go to next term
      if print_level == 0.5 or print_level >= 2:
        print("reduce gamma trace")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))
      
      reduce_two_gamma_two_partial(expression, term)
      if print_level == 0.5 or print_level >= 2:
        print("reduced two gamma two partial")
        print(term_to_LaTeX(term))
      
      TERMLEVELtrusty_reduce_eta_and_C(term)
      if print_level == 0.5 or print_level >= 2:
        print("reduce eta and C")
        print(term_to_LaTeX(term))

    gamma_chains = find_gamma_chains(term, mode = "with gamma5s") # gamma_chains list. the first sublist will correspond to tensor_eff1, tensor_eff2, etc

    tensor_effs = [] # this list will be matching the gamma_chains list in terms of tensor_effs[i] corresponds to gamma_chains[i]

    # trusting at this point in the algorithm that there are no traces left (and that gamma^2 reduced through chains)

    # now make tensor_effs
    for i in range(len(gamma_chains)):
      n = 0
      m = 0
      for tensor in gamma_chains[i]:
        if tensor.special_type == "gamma":
          n += 1
        elif tensor.special_type == "gamma5":
          m += 1
      # have n+m here because treating gamma5s like free spacetime indices
      index_structure = [["unassigned"] * (n+m),["unassigned", "unassigned"]] # does NOT matter the heights of the spacetime, as we are just using this for graphs
      tensor_eff = Gamma_Chain(index_structure, n, m)
      tensor_effs.append(tensor_eff) # not unique up to reverse ordering

    gamma_chains_set = {tensor for sublist in gamma_chains for tensor in sublist} # creating a flattened version becuase we will check it a bunch next step

    # creating tags and identifying which tensors have free indices for order
    for relation in term.relation_list:
      if relation.relation_type[:4] == "free":
        if relation.tensors[0] in gamma_chains_set: # checking for tensor in gamma_chain
          # now that we've hit, lets find which one it is in
          for i in range(len(gamma_chains)):
            if relation.tensors[0] in gamma_chains[i]:
              tensor = tensor_effs[i] # found the correct tensor_eff
              break
        else:
          tensor = relation.tensors[0]

        tag = relation.relation_type

        if tensor not in order:
          order[tensor] = [tag, tensor.special_type]
        else:
          order[tensor].append(tag)

      # contractions
      if relation.relation_type == "contraction" and not (relation.index_type == "spinorial" and relation.tensors[0].special_type in {"gamma", "gamma5"} and relation.tensors[1].special_type in {"gamma", "gamma5"}):

        gamma_chain1 = None # default
        gamma_chain2 = None

        # now we check for gamma_chain tensor_eff replacement
        if relation.tensors[0] in gamma_chains_set: # checking for tensor_eff replacement
          # now that we've hit, lets find which one it is in
          for i in range(len(gamma_chains)):
            if relation.tensors[0] in gamma_chains[i]:
              tensor1 = tensor_effs[i] # found the correct tensor_eff
              gamma_chain1 = gamma_chains[i]
              break
        else:
          tensor1 = relation.tensors[0] # could still have spacetime contractions across different gamma chains
          # this tells us nothing about tensor2 besides that it is not in the same gamma_chain as tensor1, which we excluded in the "and not" if statement

        # for tensor2 now (same code, 0 -> 1)
        if relation.tensors[1] in gamma_chains_set: # checking for tensor_eff replacement
          # now that we've hit, lets find which one it is in
          for i in range(len(gamma_chains)):
            if relation.tensors[1] in gamma_chains[i]:
              tensor2 = tensor_effs[i] # found the correct tensor_eff, know this can't be the same as tensor1 by prior trace, gamma^2 reductions
              gamma_chain2 = gamma_chains[i]
              break
        else:
          tensor2 = relation.tensors[1]

        # now that tensor_eff consideration done, let's make graphs
        graph_num = -1 # unset

        for i in range(len(graph_list)):
          if any(tensor1 in [x[0] for x in s] or tensor2 in [x[0] for x in s] for s in graph.values()): # same as s[:, 0]
            graph_num = i
            break

        if graph_num == -1:
          graph = {} # new graph
        else:
          graph = graph_list[graph_num]

        if graph.get(tensor1) is None:
          graph[tensor1] = [] # making list of edges if they don't already exist
        if graph.get(tensor2) is None:
          graph[tensor2] = []


        # setting up to/from index, with Gamma_chain consideration (to not save inner contractions)

        if gamma_chain1 is None:
          index1 = relation.indices[0]
        else:
          if relation.index_type == "spacetime":
            index1 = gamma_chain1.index(relation.tensors[0])
          else: # spinorial
            index1 = 0 if gamma_chain1.index(relation.tensors[0]) == 0 else 1
        if gamma_chain2 is None:
          index2 = relation.indices[1]
        else:
          if relation.index_type == "spacetime":
            index2 = gamma_chain2.index(relation.tensors[1])
          else: # spinorial
            index2 = 0 if gamma_chain2.index(relation.tensors[1]) == 0 else 1

        graph[tensor1].append([tensor2, relation.index_type, index1, index2]) # [..., from index, to index]
        graph[tensor2].append([tensor1, relation.index_type, index2, index1])

        if graph_num == -1:
          graph_list.append(graph)

    for tensor in order:
      order[tensor] = str(sorted(order[tensor]))

    # --------------------------------------------------------------------------------------------------------

    # merging non disjoint graphs
    length = len(graph_list)

    i = 0

    while i < length:
      length = len(graph_list)
      j = i + 1 # pairwise, order doesn't matter, don't want to hit i=j
      graph1 = graph_list[i]

      while j < length:
        graph2 = graph_list[j]
        if not graph1.keys().isdisjoint(graph2.keys()):
          # then we merge them
          for tensor in graph2:
            if graph1.get(tensor) is None:
              graph1[tensor] = graph2[tensor] # this is if the other tensor is not in here
            else:
              for edge in graph2[tensor]: # this is to complete the edges for each tensor
                if edge not in graph1[tensor]:
                  graph1[tensor].append(edge)
          del graph_list[j]
          length = len(graph_list)
        else:
          j += 1 # doesn't index because we keep checking otherwise

      i += 1

    # --------------------------------------------------------------------------------------------------------

    # now graph_list only has disjoint graphs

    for tensor in term.tensor_list: # this term_list is fine here, because doesn't deal with gamma_chains
      if tensor.special_type[:5] == "field":
        order[tensor] = tensor.special_type
        field = tensor # this is used later in build_path_level
        # should be at end of term

    tensor_eff_list = list(term.tensor_list)

    i = 0

    while i < len(tensor_eff_list):
      tensor = tensor_eff_list[i]
      if tensor in gamma_chains_set:
        tensor_eff_list.remove(tensor) # getting rid of gamma_chains to replace with tensor_eff version
        i -= 1 # deindexing because we removed one
      i += 1

    tensor_eff_list += tensor_effs # now we have our list for the graph

    # create tags
    for start_tensor in tensor_eff_list: # tensor_eff_list, not term.tensor_list
      if order.get(start_tensor) is not None:
        continue # this moves to next start_tensor

      for graph_temp in graph_list:
        if start_tensor in graph_temp.keys():
          graph = graph_temp
          break

      level = 0
      end = False
      pre_tag = [start_tensor, []]

      while end == False:

        end = build_path_level(pre_tag, graph, order, field, level)
        level += 1

      # store tensors in pre_tag, then replace with special_type later before deepsort and convert to string
      convert_and_deepsort(pre_tag)

      tag = str(pre_tag)

      # we have tags now
      order[start_tensor] = tag

    # ----------------------------------------------------------------------------------------------

    # TAGBUMPS (intergraph and outergraph)
    # tag bumps -> adding a *1 for intergaph, &1 for outergraph

    # OUTERGRAPH
    for i in range(len(tensor_eff_list)):
      bump_num = 1 # resets here because no more graphs will be the same as the previous tensor's
      tensor1 = tensor_eff_list[i]

      for j in range(i+1, len(tensor_eff_list)):
        tensor2 = tensor_eff_list[j]
        if order[tensor1] == order[tensor2]: # needs some type of bump

          graph1 = None # default
          for graph in graph_list:
            if tensor1 in graph:
              graph1 = graph # identified tensor1's graph
              break

          if graph1 is None: # tensor1, tensor2 are not in a graph
            order[tensor2] = order[tensor2] + "&" + str(bump_num)
            bump_num += 1

            # graph1 exists, so graph2 must exist for tensor2 too
          elif tensor2 not in graph1: # this is the outergraph condition
            for graph in graph_list:
              if tensor2 in graph:
                graph2 = graph # identified tensor2's graph
                break
            for tensor in graph2:
              order[tensor] = order[tensor] + "&" + str(bump_num) # bumping everything in graph 2
            bump_num += 1

    # INTERGRAPH
    for i in range(len(tensor_eff_list)):
      bump_num = 1 # resets here because no more subgraphs will be the same as the previous tensor's
      tensor1 = tensor_eff_list[i]

      # now we search the graph in progressively larger regions around it, and don't hit tensors we have already hit. if there is a tie, then it must be a 3-prong equivalence, so we tabugmp them in any order

      graph = None

      for graph_temp in graph_list:
        if tensor1 in graph_temp:
          graph = graph_temp # identified tensor1's graph
          break

      if graph == None:
        continue

      # now we use this graph to detect for integraph teagbumps
      end = False
      while end == False:

        end = True # default

        tensors_hit = [tensor1]
        tensor_to_bump = None # default

        # now we search for equivalents
        end2 = False
        while end2 == False: # go until we find an equivalent tag

          end2 = True # default
          get_out_of_loops = False

          length = len(tensors_hit)

          for j in range(length):
            tensor2 = tensors_hit[j]
            for edge in graph[tensor2]:
              if edge[0] not in tensors_hit:
                if order[edge[0]] == order[tensor1]: # needs a bump
                  # get out of all these loops and bump
                  get_out_of_loops = True
                  tensor_to_bump = edge[0]
                  break
                # else
                tensors_hit.append(edge[0])
                end2 = False
                # we bump if there were are any hits. need this have be the same order as the paths are in.
            if get_out_of_loops == True:
              break

        if tensor_to_bump is not None:
          end = False # keep searching for more

          to_bump = find_local_graph(graph1, tensor1, tensor_to_bump)

          # bump them
          for tensor in to_bump:
            order[tensor] = order[tensor] + "*" + str(bump_num) # tagbump everything in tensor2's list
          bump_num += 1

    # ----------------------------------------------------------------------------------------------

    # finals steps to get final_tensor_list
    order_arr = [list(x) for x in order.items()] # converts to lists

    sorted_order_list = sorted(order_arr, key=lambda x: x[1])

    final_tensor_eff_list = [x[0] for x in sorted_order_list] # same as sorted_order_list[:,0]

    # now we do index sorting

    # ----------------------------------------------------------------------------------------------
    # now we do the order of the indices on each tensor, the tag they are contracted to or on their free ID
    # do it tensor wise, then expand tensor_eff chains. Tagbump contracted indices to same object

    # order by lowest spacetime tag, then put the leftmost spin index on that

    space_index_list = [] # will be index sublist matching final_tensor_eff_list[i], i.e. the indices and tags on that tensor_eff
    spin_index_list = []


    for tensor in final_tensor_eff_list:
      # for space_index on this tensor
      space_index_sublist = []
      spin_index_sublist = [] # sublists to be added to it later

      gamma5_number = 0

      for space_index in range(len(tensor.index_structure[0])):
        free = True # default

        tag = None # default

        for graph in graph_list:
          if tensor in graph: # it's in here
            for edge in graph[tensor]:
              if edge[2] == space_index and edge[1] == "spacetime": # this is from_index, and relation.relation_type
                free = False
                tag = order[edge[0]]
                break # found it
            break # found the tag or found that it's not here

        if free == True: # it's a free index
          if tensor.special_type[:11] == "gamma_chain": # checking for tensor in gamma_chain
            gamma_chain = gamma_chains[tensor_effs.index(tensor)]
            for relation in term.relation_list:
              if relation.relation_type[:4] == "free" and relation.index_type == "spacetime": # filters
                if relation.tensors[0] is gamma_chain[space_index]:
                  tag = relation.relation_type
                  break
          else:
            for relation in term.relation_list:
              if relation.relation_type[:4] == "free" and relation.index_type == "spacetime": # filters
                if relation.tensors[0] is tensor and relation.indices[0] == space_index:
                  tag = relation.relation_type   # there is only one tensor in relation.tensors for a free index relation
                  break # found it

        if tag is None and tensor.special_type[:11] == "gamma_chain":
          tag = "gamma5 " + str(gamma5_number)
          gamma5_number += 1

        space_index_sublist.append([space_index, tag])
        # tag of that index = tag of the contracted object
        # (if its a gamma_chain, go the ith gamma to get the relevant spacetime index or the edge ones if its a spinorial)
      for spin_index in range(len(tensor.index_structure[1])):
        free = True # default

        for graph in graph_list:
          if tensor in graph: # it's in here
            for edge in graph[tensor]:
              if edge[2] == spin_index and edge[1] == "spinorial": # this is from_index, and relation.relation_type
                free = False
                tag = order[edge[0]]
                break # found it
            break # found the tag or found that it's not here

        if free == True: # its a free index
          if tensor.special_type[:11] == "gamma_chain": # checking for tensor in gamma_chain
            gamma_chain = gamma_chains[tensor_effs.index(tensor)]
            temp_index = 0 if spin_index == 0 else -1 # it will either be the first or the last of the chain
            for relation in term.relation_list:
              if relation.relation_type[:4] == "free" and relation.index_type == "spinorial": # filters
                if relation.tensors[0] is gamma_chain[temp_index]:
                  tag = relation.relation_type
                  break
          else:
            for relation in term.relation_list:
              if relation.relation_type[:4] == "free" and relation.index_type == "spinorial": # filters
                if relation.tensors[0] is tensor and relation.indices[0] == spin_index:
                  tag = relation.relation_type   # there is only one tensor in relation.tensors for a free index relation
                  break # found it

        # not same as space, since only corners

        spin_index_sublist.append([spin_index, tag])

      # can't order them yet, have to tagbump first
      space_index_list.append(space_index_sublist)
      spin_index_list.append(spin_index_sublist)

    # ----------------------------------------------------------------------------------------------

    # tagbump the indices!
    # have to tagbump the same contraction

    # tagbump part.
    # only contracted indices will to tagbump, so can use graph. only have to tagbump across two tensors, and for space and spin seperately
    # index tagbump is #1. << 1 is bump_num here

    # INDEX tagbumps (spin)
    for tensor_num in range(len(spin_index_list)):
      spin_index_sublist = spin_index_list[tensor_num]

      for i in range(len(spin_index_sublist)):
        equivalent_tags = []
        tag1 = spin_index_sublist[i][1] # comparing if two tags are the same
        for j in range(i+1, len(space_index_sublist)):  # don't want to tagbump self
          tag2 = spin_index_sublist[j][1]
          if tag1 == tag2:
            equivalent_tags.append(spin_index_sublist[j]) # tagbump the second element of these sublists (list is mutable)

        # tagbumps
        tensor = final_tensor_eff_list[tensor_num] # corresponding tensor
        other_contraction_bump_num = 0
        self_contraction_bump_num = 0
        temp = 0

        # we go through each, and we bump the tag and its contracted index on the other tensor
        while temp < len(equivalent_tags):
          equivalent_tag_pair = equivalent_tags[temp]
          temp += 1

          for graph in graph_list:
            if tensor in graph: # it's in here
              for edge in graph[tensor]:
                if edge[2] == equivalent_tag_pair[0] and edge[1] == "spinorial": # this is from_index, and relation.relation_type
                  tensor2 = edge[0] # tensor with the other index we need to tagbump (contracted)

                  if tensor is tensor2: # then this is self contraction

                    equivalent_tag_pair[1] = equivalent_tag_pair[1] + "!" + str(self_contraction_bump_num)
                    spin_index_list[tensor_num][edge[3]][1] = spin_index_list[tensor_num][edge[3]][1] + "?" + str(self_contraction_bump_num)
                    self_contraction_bump_num += 1

                    # deleting the other self contracted index:
                    for k in range(len(equivalent_tags)):
                      if equivalent_tags[k][0] == edge[3]: # this checks to find the index_num of the other self contracted index
                        del equivalent_tags[k]
                        break

                  else:
                    tensor2_num = final_tensor_eff_list.index(tensor2) # this is where the tensor2 is in final_tensor_eff_list

                    # now we have to find the corresponding index, edge[3] is to_index
                    equivalent_tag_pair[1] = equivalent_tag_pair[1] + "#" + str(other_contraction_bump_num)
                    spin_index_list[tensor2_num][edge[3]][1] = spin_index_list[tensor2_num][edge[3]][1] + "#" + str(other_contraction_bump_num)
                    other_contraction_bump_num += 1
                  # this tagbump done!
                  break # found it
              break # found it

    # INDEX tagbumps (space)
    for tensor_num in range(len(space_index_list)):
      space_index_sublist = space_index_list[tensor_num]
      for i in range(len(space_index_sublist)):
        equivalent_tags = []
        tag1 = space_index_sublist[i][1] # comparing if two tags are the same
        for j in range(i+1, len(space_index_sublist)): # don't want to tagbump self
          tag2 = space_index_sublist[j][1]
          if tag1 == tag2:
            equivalent_tags.append(space_index_sublist[j]) # tagbump the second element of these sublists (list is mutable)

        # tagbumps
        tensor = final_tensor_eff_list[tensor_num] # corresponding tensor
        other_contraction_bump_num = 0
        self_contraction_bump_num = 0
        temp = 0

        # we go through each, and we bump the tag and its contracted index on the other tensor
        while temp < len(equivalent_tags):
          equivalent_tag_pair = equivalent_tags[temp]
          temp += 1

          for graph in graph_list:
            if tensor in graph: # it's in here
              for edge in graph[tensor]:
                if edge[2] == equivalent_tag_pair[0] and edge[1] == "spacetime": # this is from_index, and relation.relation_type
                  tensor2 = edge[0] # tensor with the other index we need to tagbump (contracted)

                  if tensor is tensor2: # then this is self contraction

                    equivalent_tag_pair[1] = equivalent_tag_pair[1] + "!" + str(self_contraction_bump_num)
                    space_index_list[tensor_num][edge[3]][1] = space_index_list[tensor_num][edge[3]][1] + "?" + str(self_contraction_bump_num)
                    self_contraction_bump_num += 1

                    # deleting the other self contracted index:
                    for k in range(len(equivalent_tags)):
                      if equivalent_tags[k][0] == edge[3]: # this checks to find the index_num of the other self contracted index
                        del equivalent_tags[k]
                        break

                  else:
                    tensor2_num = final_tensor_eff_list.index(tensor2) # this is where the tensor2 is in final_tensor_eff_list

                    # now we have to find the corresponding index, edge[3] is to_index
                    equivalent_tag_pair[1] = equivalent_tag_pair[1] + "#" + str(other_contraction_bump_num)
                    space_index_list[tensor2_num][edge[3]][1] = space_index_list[tensor2_num][edge[3]][1] + "#" + str(other_contraction_bump_num)

                    other_contraction_bump_num += 1
                  # this tagbump done!
                  break # found it
              break # found it


    # ----------------------------------------------------------------------------------------------

    space_index_list2 = [sorted(sublist, key=lambda x: x[1]) for sublist in space_index_list]
    spin_index_list2 = [sorted(sublist, key=lambda x: x[1]) for sublist in spin_index_list]
    # now we have a way to reorder the tensor indices

    flip_list = [None]*len(final_tensor_eff_list)

    coeff_factor = 1 # this is held for a while

    for i in range(len(final_tensor_eff_list)):
      tensor = final_tensor_eff_list[i]

      # spinorial first
      temp_list = [x[0] for x in spin_index_list[i]]
      final_list = [x[0] for x in spin_index_list2[i]]

      for switching_index in range(len(temp_list)): # switching_index is the one we are going to move into position, e.g. 1 2 3 or 4. we search for it now
        index_pos1 = temp_list.index(switching_index)
        index_pos2 = final_list.index(switching_index)
        if index_pos1 != index_pos2: # then it is NOT already in position

          if tensor.special_type[:11] == "gamma_chain":
            flip_list[i] = True
            # flipping order of chain or not

          else: # not on a tensor_eff
            sym = symmetry(tensor, index_pos1, index_pos2, "spinorial")
            if sym is not None:
              temp = 1 if sym == "sym" else -1
              coeff_factor *= temp

              # reassiging indices

              index1 = temp_list[index_pos1]
              index2 = temp_list[index_pos2]
              temp_list[index_pos1] = index2
              temp_list[index_pos2] = index1

              # reassiging heights

              index1_height = tensor.index_structure[1][index_pos1]
              index2_height = tensor.index_structure[1][index_pos2]
              tensor.index_structure[1][index_pos1] = index2_height
              tensor.index_structure[1][index_pos2] = index1_height

              # reassiging relations

              for relation in term.relation_list:
                if relation.index_type == "spinorial":
                  for l in range(len(relation.tensors)):
                    if relation.tensors[l] is tensor and relation.indices[l] == index_pos1:
                      relation.indices[l] = index_pos2
                      # if it's a contraction between j and k do nothing, or we can switch it does not matter
                      # only problem is we don't want to hit j's we turned into k's and turn them back into j's, by only going over each index once it's fine
                    elif relation.tensors[l] is tensor and relation.indices[l] == index_pos2:
                      relation.indices[l] = index_pos1

      # spacetime second

      temp_list = [x[0] for x in space_index_list[i]]
      final_list = [x[0] for x in space_index_list2[i]]

      for switching_index in range(len(temp_list)): # switching_index is the one we are gonna move into position, e.g. 1 2 3 or 4. we search for it now
        index_pos1 = temp_list.index(switching_index)
        index_pos2 = final_list.index(switching_index)
        if index_pos1 != index_pos2: # then it is NOT already in position
          if tensor.special_type[:11] != "gamma_chain":
            sym = symmetry(tensor, index_pos1, index_pos2, "spacetime")
            if sym is not None:
              temp = 1 if sym == "sym" else -1
              coeff_factor *= temp

              # rearranging indices

              index1 = temp_list[index_pos1]
              index2 = temp_list[index_pos2]
              temp_list[index_pos1] = index2
              temp_list[index_pos2] = index1

              # reassiging heights

              index1_height = tensor.index_structure[0][index_pos1]
              index2_height = tensor.index_structure[0][index_pos2]
              tensor.index_structure[0][index_pos1] = index2_height
              tensor.index_structure[0][index_pos2] = index1_height

              # reassiging relations

              for relation in term.relation_list:
                if relation.index_type == "spacetime":
                  for l in range(len(relation.tensors)):
                    if relation.tensors[l] is tensor and relation.indices[l] == index_pos1:
                      relation.indices[l] = index_pos2
                      # if it's a contraction between j and k do nothing, or we can switch it does not matter
                      # only problem is we don't want to hit j's we turned into k's and turn them back into j's, by only going over each index once its fine
                    elif relation.tensors[l] is tensor and relation.indices[l] == index_pos2:
                      relation.indices[l] = index_pos1

    term.coefficient *= coeff_factor # want clifford terms to have the correct sign since we did switches

    # apply Clifford for spacetime

    for i in range(len(final_tensor_eff_list)):
      tensor = final_tensor_eff_list[i]

      if tensor.special_type[:11] == "gamma_chain":
        gamma_chain = gamma_chains[tensor_effs.index(tensor)]

        from_list = [x[0] for x in space_index_list[i]] # same as space_index_list[i][:,0]
        if flip_list[i] == True:
          to_list = [x[0] for x in space_index_list2[i]][::-1] # same as space_index_list2[i][:,0][::-1]
        else:
          to_list = [x[0] for x in space_index_list2[i]] # same as space_index_list2[i][:,0]

        # now run routine swap until from_list matches to list
        # take the first one and flip until it reaches desired spot, then second, etc

        for target_index in range(len(to_list)):
          index_pos1 = from_list.index(target_index)
          index_pos2 = to_list.index(target_index)
          if index_pos1 < index_pos2:
            for swap_index in range(index_pos1, index_pos2):
              gammaL = gamma_chain[swap_index]
              gammaR = gamma_chain[swap_index+1]

              for relation in term.relation_list:
                if relation.relation_type == "contraction" and relation.index_type == "spinorial": # know gamma^2 won't exist anymore
                  if (relation.tensors[0] is gammaL and relation.tensors[1] is gammaR) or (relation.tensors[1] is gammaL and relation.tensors[0] is gammaR): # trusting no space contraction between these as well
                    relationA = relation
                    break # found it

              apply_Clifford(expression, term, gammaL, gammaR, relationA)

              gamma_chain[swap_index] = term.tensor_list[-2]
              gamma_chain[swap_index+1] = term.tensor_list[-1] # the last appended tensor will be new tensor to put here, from format of how apply_Clifford works

              # now switching in the from_list
              index1 = from_list[swap_index]
              index2 = from_list[swap_index+1]
              from_list[swap_index] = index2
              from_list[swap_index+1] = index1

          elif index_pos1 > index_pos2:
            for swap_index in range(index_pos1, index_pos2, -1):
              gammaL = gamma_chain[swap_index-1]
              gammaR = gamma_chain[swap_index]

              for relation in term.relation_list:
                if relation.relation_type == "contraction" and relation.index_type == "spinorial":
                  if (relation.tensors[0] is gammaL and relation.tensors[1] is gammaR) or (relation.tensors[1] is gammaL and relation.tensors[0] is gammaR): # trusting no space contraction between these as well
                    relationA = relation
                    break # found it

              apply_Clifford(expression, term, gammaL, gammaR, relationA)

              gamma_chain[swap_index] = term.tensor_list[-1]
              gamma_chain[swap_index-1] = term.tensor_list[-2] # the second to last appended tensor will be here in the chain

              # now switching in the from_list
              index1 = from_list[swap_index-1]
              index2 = from_list[swap_index]
              from_list[swap_index-1] = index2
              from_list[swap_index] = index1

    # now we fix spinorial convention
    # creating dummy spinorial contraction list

    # also details gamma5 swtiching
    gamma5s_to_swap = set()

    coeff_factor = 1

    spin_contractions = []
    for relation in term.relation_list:
      if relation.relation_type == "contraction" and relation.index_type == "spinorial":
        spin_contractions.append(relation)

    for tensor in final_tensor_eff_list:
      if tensor.special_type[:11] == "gamma_chain":
        gamma_chain = gamma_chains[tensor_effs.index(tensor)]

        flip = flip_list[final_tensor_eff_list.index(tensor)]
        if flip == True:
          gamma_chain = gamma_chain[::-1]

        for gamma in gamma_chain: # doing it in order
          for spin_index in range(len(gamma.index_structure[1])):
            for i in range(len(spin_contractions)):
              spin_contraction = spin_contractions[i]
              leave_loop = False
              for j in range(len(spin_contraction.tensors)):
                if spin_contraction.tensors[j] is gamma and spin_contraction.indices[j] == spin_index and spin_contraction.tensors[1-j] in gamma_chain:

                  tensor2 = spin_contraction.tensors[1-j]
                  index2 = spin_contraction.indices[1-j]

                  if gamma.index_structure[1][spin_index] == "down":
                    gamma.index_structure[1][spin_index] = "up"
                    # now switching other one
                    tensor2.index_structure[1][index2] = "down"
                    # flipping index height does this:
                    coeff_factor *= -1
                    # if it is down we do nothing

                  # now gamma5 switching if necessary
                  if gamma.special_type == "gamma5":
                    if spin_index == 0:
                      gamma5s_to_swap.add(gamma)

                  if tensor2.special_type == "gamma5":
                    if index2 == 1: # this is not convention, we want to swap indices
                      gamma5s_to_swap.add(tensor2)

                  del spin_contractions[i] # remove it since we already fixed the convention for this pair
                  leave_loop = True
                  break
              if leave_loop == True:
                break # go to next one
      else:
        for spin_index in range(len(tensor.index_structure[1])):
          for i in range(len(spin_contractions)):
            spin_contraction = spin_contractions[i]
            leave_loop = False
            for j in range(len(spin_contraction.tensors)):
              if spin_contraction.tensors[j] is tensor and spin_contraction.indices[j] == spin_index:
                if tensor.index_structure[1][spin_index] == "up":
                  tensor.index_structure[1][spin_index] = "down"
                  # now switching other one
                  tensor2 = spin_contraction.tensors[1-j]
                  index2 = spin_contraction.indices[1-j]
                  tensor2.index_structure[1][index2] = "up"
                  # flipping index height does this:
                  coeff_factor *= -1
                  # if it is down we do nothing
                del spin_contractions[i] # remove it since we already fixed the convention for this pair
                leave_loop = True
                break
            if leave_loop == True:
              break # go to next one

    init_spin_dict = {}
    for gamma5 in gamma5s_to_swap:
      init_spin_dict[gamma5] = gamma5.index_structure[1].copy() # if we don't copy, will get errors

    for relation in term.relation_list:
      for i in range(len(relation.tensors)):
        tensor = relation.tensors[i]
        index = relation.indices[i]
        if tensor in gamma5s_to_swap:
          relation.indices[i] = 1 - index # swaps 1 and 0
          tensor.index_structure[1][1 - index] = init_spin_dict[tensor][index]


    coeff_factor *= (-1)**len(gamma5s_to_swap) # 1 swap is 1 sign flip, 2 is 2, etc


    term.coefficient *= coeff_factor # multiplying

    tag_sublist = []

    for x in sorted_order_list:
      tag_sublist.append(x[1])

    tag_list.append(tag_sublist)

    for i in range(len(term.tensor_list)):
      tensor = term.tensor_list[i]
      if tensor.special_type == "partial":
        term.tensor_list.append(term.tensor_list.pop(i))

    for i in range(len(term.tensor_list)):
      tensor = term.tensor_list[i]
      if tensor.special_type[:5] == "field":
        term.tensor_list.append(term.tensor_list.pop(i))

  if print_level == 0.5 or print_level >= 1.5:
    print("expression post tagmaker:")
    print("0 = " + expression_to_LaTeX(expression))

  return tag_list # outputting this for matcher

# trusty that expand_gamma5s has already been applied

def trusty_reduce_long_chains(expression, term):

  end = False
  while end == False:

    end = True # default

    gamma_chains = find_gamma_chains(term)

    gamma_chain = None

    for gamma_chain_temp in gamma_chains:
      if len(gamma_chain_temp) >= 4:
        end = False
        gamma_chain = gamma_chain_temp
        break

    if gamma_chain is None:
      break # redundant but makes code clearer to read

    # now we apply the identity to the gamma chain

    reduce_chain_once(gamma_chain, expression, term)


def reduce_chain_once(gamma_chain, expression, term):

  gamma1 = gamma_chain[0]
  gamma2 = gamma_chain[1]
  gamma3 = gamma_chain[2]
  # want output to give us the one with one gamma as the first term

  old_tensors = {gamma,gamma2,gamma3}

  old_relations = set()


  for relation in term.relation_list:
    if relation.relation_type == "contraction" and relation.index_type == "spinorial":
      for i in range(len(relation.tensors)):
        tensor1 = relation.tensors[i]
        tensor2 = relation.tensors[1-i]
        if tensor1 is gamma1 and tensor2 is gamma2:
          old_relations.add(relation)
          a_index = relation.indices[i]
          break
        if tensor1 is gamma2 and tensor2 is gamma3:
          old_relations.add(relation)
          b_index = index.indices[1-i]
          break

  # now creating the rest of the items for TERMLEVELapply_identity

  mu = gamma1.index_structure[0][0]
  nu = gamma2.index_structure[0][0]
  rho = gamma3.index_structure[0][0]
  a = gamma1.index_structure[1][a_index]
  b = gamma3.index_structure[1][b_index]

  new_tensors_list = []
  new_relations_list = []
  tensor_index_dict = {}

  # term 1
  eta = Eta([[mu,nu],[]])
  gamma = Gamma_mu([[rho],[a,b]])
  new_tensors = [eta, gamma]

  new_relations = [] # will be used for first 3 terms

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(gamma1, 0, 0, "spacetime")] = [eta, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(gamma2, 0, 0, "spacetime")] = [eta, 1]
  tensor_index_dict[(gamma3, 0, 0, "spacetime")] = [gamma, 0]

  tensor_index_dict[(gamma1, a_index, 0, "spinorial")] = [gamma, 0]
  tensor_index_dict[(gamma3, b_index, 0, "spinorial")] = [gamma, 1]

  # term 2
  eta = Eta([[mu,rho],[]])
  gamma = Gamma_mu([[nu],[a,b]])
  new_tensors = [eta, gamma]

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(gamma1, 0, 1, "spacetime")] = [eta, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(gamma2, 0, 1, "spacetime")] = [gamma, 0]
  tensor_index_dict[(gamma3, 0, 1, "spacetime")] = [eta, 1]

  tensor_index_dict[(gamma1, a_index, 1, "spinorial")] = [gamma, 0]
  tensor_index_dict[(gamma3, b_index, 1, "spinorial")] = [gamma, 1]

  # term 2
  eta = Eta([[nu,rho],[]])
  gamma = Gamma_mu([[mu],[a,b]])
  new_tensors = [eta, gamma]

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(gamma1, 0, 2, "spacetime")] = [gamma, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(gamma2, 0, 2, "spacetime")] = [eta, 0]
  tensor_index_dict[(gamma3, 0, 2, "spacetime")] = [eta, 1]

  tensor_index_dict[(gamma1, a_index, 2, "spinorial")] = [gamma, 0]
  tensor_index_dict[(gamma3, b_index, 2, "spinorial")] = [gamma, 1]

  # term3
  Levi = Levi_Civita([[mu,nu,rho,"up"],[]])
  gamma = Gamma_mu([["down"],[a,"up"]])
  gamma5 = Gamma5([[],["down", b]])
  new_tensors = [Levi, gamma, gamma5]

  relation1 = Contraction([Levi, gamma], [3, 0], "spacetime")
  relation2 = Contraction([gamma, gamma5], [1, 0], "spinorial")
  new_relations = [relation1, relation2]

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(gamma1, 0, 3, "spacetime")] = [Levi, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(gamma2, 0, 3, "spacetime")] = [Levi, 1]
  tensor_index_dict[(gamma3, 0, 3, "spacetime")] = [Levi, 2]

  tensor_index_dict[(gamma1, a_index, 3, "spinorial")] = [gamma, 0]
  tensor_index_dict[(gamma3, b_index, 3, "spinorial")] = [gamma5, 1]


  # now doing coeff_factor_list
  overall_factor = 1 # spinorial height convention
  if gamma1.index_structure[1][1 - a_index] == "down":
    overall_factor *= -1
  if gamma3.index_structure[1][1 - b_index] == "up":
    overall_factor *= -1

  coeff_factor_list = [overall_factor, -1 * overall_factor, overall_factor, -1j * overall_factor]


  output_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
  expression.term_list = expression.term_list + output_terms


def trusty_reduce_eta_and_C(expression):
  i = 0
  while i < len(expression.term_list):
    term = expression.term_list[i]
    i += 1
    zeroed = TERMLEVELtrusty_reduce_eta_and_C(term) # remove term
    if zeroed == True:
      expression.term_list.remove(term)
      i -= 1 # deindexing because we removed one




def TERMLEVELtrusty_reduce_eta_and_C(term):#, expression):

  # let's identify eta and C's that can be reduced

  num = 0
  while num < len(term.tensor_list):

    num += 1

    found_one = False

    for temp in range(len(term.relation_list)):
      relation = term.relation_list[temp]
      if relation.relation_type == "contraction":
        for i in range(len(relation.tensors)):
          tensor = relation.tensors[i]
          if tensor.special_type == "metric":

            found_one = True
            eta = tensor

            # set up objects for TERMLEVELapply_identity
            tensor2 = relation.tensors[1-i] # the other tensor

            old_tensors = {tensor2, eta}
            old_relations = {relation}

            if relation.tensors[0] is relation.tensors[1]: # self contraction
              new_tensors_list = [[]]
              new_relations_list = [[]]
              tensor_index_dict = {}
              coeff_factor_list = [4]
              break # won't continue the code


            index2 = relation.indices[1-i]
            to_index = 1 - relation.indices[i] # the other index on eta

            tensor2_copy = tensor2.copy()
            to_height = eta.index_structure[0][to_index]# the other index on eta
            tensor2_copy.index_structure[0][index2] = to_height

            new_tensors_list = [[tensor2_copy]]
            new_relations_list = [[]]

            tensor_index_dict = {}
            tensor_index_dict[(eta, to_index, 0, "spacetime")] = [tensor2_copy, index2] # tuple is (tensor, index, term_num, relation.relation_type)
            for j in range(len(tensor2_copy.index_structure[0])):
              if j != index2:
                tensor_index_dict[(tensor2, j, 0, "spacetime")] = [tensor2_copy, j]
            for j in range(len(tensor2_copy.index_structure[1])):
              tensor_index_dict[(tensor2, j, 0, "spinorial")] = [tensor2_copy, j]

            coeff_factor_list = [1]

            break # found an eta to reduce, stop searching

          elif tensor.special_type == "C_matrix":

            found_one = True
            C = tensor

            # set up objects for TERMLEVELapply_identity
            tensor2 = relation.tensors[1-i] # the other tensor

            old_tensors = {tensor2, C}
            old_relations = {relation}

            if relation.tensors[0] is relation.tensors[1]: # self contraction
              return True # this term is zeroed

            index2 = relation.indices[1-i]
            to_index = 1 - relation.indices[i] # the other index on eta

            tensor2_copy = tensor2.copy()
            to_height = C.index_structure[1][to_index]# the other index on eta
            tensor2_copy.index_structure[1][index2] = to_height

            new_tensors_list = [[tensor2_copy]]
            new_relations_list = [[]]

            tensor_index_dict = {}
            tensor_index_dict[(C, to_index, 0, "spinorial")] = [tensor2_copy, index2] # tuple is (tensor, index, term_num, relation.relation_type)
            for j in range(len(tensor2_copy.index_structure[0])):
              tensor_index_dict[(tensor2, j, 0, "spacetime")] = [tensor2_copy, j]
            for j in range(len(tensor2_copy.index_structure[1])):
              if j != index2:
                tensor_index_dict[(tensor2, j, 0, "spinorial")] = [tensor2_copy, j]


            # with C there is an up/down sign convention!
            C_index1 = relation.indices[i]
            if (C.index_structure[1][C_index1] == "down" and C_index1 == 0) or (C.index_structure[1][C_index1] == "up" and C_index1 == 1):
              sign = 1
            else:
              sign = -1

            coeff_factor_list = [sign]

            break # found a C to reduce, stop searching

        if found_one == True:
          break # leave this loop because we want to apply_identity now

    if found_one == True:
      TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
      # will have no output_terms
      num -= 1 # deindex since we just removed a tensor

  return False # zeroed = False

# trusty because assumes no excess etas and no symmetrizers

def trusty_use_Fierz(expression):

  trusty_reduce_eta_and_C(expression)

  # now we do the gammas Fierz
  term_num = 0
  num_of_terms = len(expression.term_list)
  while term_num < num_of_terms: # don't want to double hit terms with the same current_ID

    term = expression.term_list[term_num]
    term_num += 1

    hits = 0
    for tensor in term.tensor_list:
      if tensor.special_type == "partial":
        hits += 1
        break

    if hits != 1:
      continue # won't activate for things with multiple partials in a term

    gamma_chains = find_gamma_chains(term, mode = "with gamma5s", singletons = True)

    gamma_chain = None

    # identifying the gamma with mu (no sigma so garuanteed to exist)
    for relation in term.relation_list:
      if relation.index_type == "spacetime" and relation.relation_type == "contraction":
        for gamma_chain_temp in gamma_chains:
          if relation.tensors[0].special_type == "partial" and relation.tensors[1] in gamma_chain_temp:
            gamma_chain = gamma_chain_temp
            break # found it!
          elif relation.tensors[1].special_type == "partial" and relation.tensors[0] in gamma_chain_temp:
            gamma_chain = gamma_chain_temp
            break # found it!
        if gamma_chain is not None:
          break

    if gamma_chain is None:
      continue

    # if gamma chain has the a we put the b on it. if gamma chain b we put the a on it. if niether continue

    for relation in term.relation_list:
      if relation.relation_type == "free1":
        tensor_a = relation.tensors[0]
        index_a = relation.indices[0]
      if relation.relation_type == "free2":
        tensor_b = relation.tensors[0]
        index_b = relation.indices[0]

    gamma_chain_set = set(gamma_chain)

    tensor1 = tensor_a
    index1 = index_a
    tensor2 = tensor_b
    index2 = index_b

    output_terms = TERMLEVELtrusty_use_Fierz(term, tensor1, index1, tensor2, index2)
    expression.term_list = expression.term_list + output_terms

def TERMLEVELtrusty_use_Fierz(term, tensor1, index1, tensor2, index2):

  # assuming we already found tensor1, tensor2

  # building up objects for TERMLEVELapply_identity

  old_tensors = {tensor1, tensor2}
  old_relations = set()

  # this is a lot of terms since it's a Fierz
  new_tensors_list = []
  new_relations_list = []
  tensor_index_dict = {}

  # gamma term
  a = tensor1.index_structure[1][index1] # going to switch a and c
  b = tensor2.index_structure[1][index2]
  c = "up" if a == "down" else "down"
  d = "up" if b == "down" else "down"

  gammaL = Gamma_mu([["up"],[c,d]]) # do I need to be careful about spawning this C? ###############
  gammaR = Gamma_mu([["down"],[a,b]])
  tensor1_copy = tensor1.copy() # heights preserved
  tensor2_copy = tensor2.copy()
  new_tensors = [tensor1_copy, tensor2_copy, gammaL, gammaR]

  relation1 = Contraction([tensor1_copy, gammaL], [index1, 0], "spinorial")
  relation2 = Contraction([tensor2_copy, gammaL], [index2, 1], "spinorial")
  relation3 = Contraction([gammaL, gammaR], [0, 0], "spacetime")
  new_relations = [relation1, relation2, relation3]

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(tensor1, index1, 0, "spinorial")] = [gammaR, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(tensor2, index2, 0, "spinorial")] = [gammaR, 1] # first one was a, this is c

  for i in range(len(tensor1_copy.index_structure[0])): # reassigning relations to the new tensor copy
    tensor_index_dict[(tensor1, i, 0, "spacetime")] = [tensor1_copy, i]
  for i in range(len(tensor1_copy.index_structure[1])):
    if i != index1: # index1 shifted off tensor1, don't reattach it
      tensor_index_dict[(tensor1, i, 0, "spinorial")] = [tensor1_copy, i]
  for i in range(len(tensor2_copy.index_structure[0])):
    tensor_index_dict[(tensor2, i, 0, "spacetime")] = [tensor2_copy, i]
  for i in range(len(tensor2_copy.index_structure[1])):
    if i != index2:
      tensor_index_dict[(tensor2, i, 0, "spinorial")] = [tensor2_copy, i]

  # sigma term
  a = tensor1.index_structure[1][index1] # going to switch a and c
  b = tensor2.index_structure[1][index2]
  c = "up" if a == "down" else "down"
  d = "up" if b == "down" else "down"

  sigmaL = Sigma([["up","up"],[c,d]]) # do I need to be careful about spawning this C? ###############
  sigmaR = Sigma([["down","down"],[a,b]])
  tensor1_copy = tensor1.copy() # heights preserved
  tensor2_copy = tensor2.copy()
  new_tensors = [tensor1_copy, tensor2_copy, sigmaL, sigmaR]

  relation1 = Contraction([tensor1_copy, sigmaL], [index1, 0], "spinorial")
  relation2 = Contraction([tensor2_copy, sigmaL], [index2, 1], "spinorial")
  relation3 = Contraction([sigmaL, sigmaR], [0, 0], "spacetime")
  relation4 = Contraction([sigmaL, sigmaR], [1, 1], "spacetime")
  new_relations = [relation1, relation2, relation3, relation4]

  new_tensors_list.append(new_tensors)
  new_relations_list.append(new_relations)

  tensor_index_dict[(tensor1, index1, 1, "spinorial")] = [sigmaR, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(tensor2, index2, 1, "spinorial")] = [sigmaR, 1] # first one was a, this is c

  for i in range(len(tensor1_copy.index_structure[0])): # reassigning relations to the new tensor copy
    tensor_index_dict[(tensor1, i, 1, "spacetime")] = [tensor1_copy, i]
  for i in range(len(tensor1_copy.index_structure[1])):
    if i != index1: # index1 shifted off tensor1, don't reattach it
      tensor_index_dict[(tensor1, i, 1, "spinorial")] = [tensor1_copy, i]
  for i in range(len(tensor2_copy.index_structure[0])):
    tensor_index_dict[(tensor2, i, 1, "spacetime")] = [tensor2_copy, i]
  for i in range(len(tensor2_copy.index_structure[1])):
    if i != index2:
      tensor_index_dict[(tensor2, i, 1, "spinorial")] = [tensor2_copy, i]

  # now we do coeff_factor_list (convention of factor = 1)

  coeff_factor_list = [-1/2, -1/4] # coefficients after doing (ab)

  output_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

  return output_terms # these will be appended to expression.term_list

def trusty_reduce_gamma_squared(expression, term):

  end = False

  while end == False:
    end = True # default

    TERMLEVELtrusty_reduce_eta_and_C(term)

    gamma_chains = find_gamma_chains(term, mode = "with gamma5s")

    if len(gamma_chains) == 0:
      return # no chance to reduce

    gamma_square_exists = False # default, for if there are no gamma^2 in this term

    for gamma_chain_temp in gamma_chains:
      tensor1 = None
      for tensor in gamma_chain_temp:
        if tensor.special_type == "gamma5" and tensor1 is None:
          tensor1 = tensor
        elif tensor.special_type == "gamma5" and tensor1 is not None:
          tensor2 = tensor
          gamma_chain = gamma_chain_temp
          gamma_square_exists = True
          break # gamma^5 exists
      if gamma_square_exists == True:
        break

    if gamma_square_exists == False:
      # identifying the gamma_chain with
      for relation in term.relation_list:
        if relation.index_type == "spacetime" and relation.relation_type == "contraction" and relation.tensors[0].special_type == "gamma" and relation.tensors[1].special_type == "gamma":
          tensor1_temp = relation.tensors[0]
          tensor2_temp = relation.tensors[1]

          # checking if this is a gamma^2 contraction
          for gamma_chain_temp in gamma_chains:
            if tensor1_temp in gamma_chain_temp and tensor2_temp in gamma_chain_temp:
              gamma_square_exists = True
              gamma_chain = gamma_chain_temp
              tensor1 = tensor1_temp
              tensor2 = tensor2_temp
              space_relation = relation
              break # found it

          if gamma_square_exists == True:
            break

    if gamma_square_exists == False:
      return # no hope for this term

    num1 = -1 # unset
    num2 = -1 # unset

    for i in range(len(gamma_chain)):
      if tensor1 is gamma_chain[i]:
        num1 = i
      elif tensor2 is gamma_chain[i]:
        num2 = i

    if abs(num1 - num2) == 1:
      # now we apply gamma^2 identity since they are next to each other

      if num1 < num2:
        gammaL = gamma_chain[num1]
        gammaR = gamma_chain[num2]
      else:
        gammaL = gamma_chain[num2]
        gammaR = gamma_chain[num1]

      old_tensors = {gammaL, gammaR}

      for relationA in term.relation_list:
        if relationA.index_type == "spinorial" and relationA.relation_type == "contraction":
          if relationA.tensors[0] is gammaL and  relationA.tensors[1] is gammaR:
            a = gammaL.index_structure[1][1 - relationA.indices[0]]
            c1 = gammaL.index_structure[1][relationA.indices[0]]     # for index_structure and change_height_factor
            c2 = gammaR.index_structure[1][relationA.indices[1]]
            b = gammaR.index_structure[1][1 - relationA.indices[1]]
            num_a = 1 - relationA.indices[0]     # used for tensor_index_dict
            num_b = 1 - relationA.indices[1]
            spin_relation = relationA
            break
          elif relationA.tensors[1] is gammaL and relationA.tensors[0] is gammaR:
            a = gammaL.index_structure[1][1 - relationA.indices[1]]
            c1 = gammaL.index_structure[1][relationA.indices[1]]
            c2 = gammaR.index_structure[1][relationA.indices[0]]
            b = gammaR.index_structure[1][1 - relationA.indices[0]]
            num_a = 1 - relationA.indices[1]     # used for tensor_index_dict
            num_b = 1 - relationA.indices[0]
            spin_relation = relationA
            break

      # change_height_factor dependent on contracted index of gammaL
      if c1 == "up":
        change_height_factor = 1
      else:
        change_height_factor = -1

      # building up args for TERMLEVELapply_identity
      if gammaL.special_type == "gamma5":
        old_relations = {spin_relation}
        coeff_factor = 1
        if num_a == 1:
          coeff_factor *= -1 # price of flipping gamma5
        if num_b == 0:
          coeff_factor *= -1 # price of flipping gamma5
      else:
        old_relations = {space_relation, spin_relation}
        coeff_factor = 4

      C = C_matrix([[],[a,b]])

      new_tensors_list = [[C]]
      new_relations_list = [[]] # no new ones

      tensor_index_dict = {}
      tensor_index_dict[(gammaL, num_a, 0, "spinorial")] = [C, 0] # tuple is (tensor, index, term_num, relation.relation_type)
      tensor_index_dict[(gammaR, num_b, 0, "spinorial")] = [C, 1]

      coeff_factor_list = [coeff_factor * change_height_factor]

      TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

      end = False
    elif num1 < num2: # then we've found a gamma^2
      gammaL = gamma_chain[num2-1]
      gammaR = gamma_chain[num2]
      end = False
    elif num2 < num1:
      gammaL = gamma_chain[num1-1]
      gammaR = gamma_chain[num1]
      end = False
      # does nothing if num1 == num2 (i.e. if still unset)
    # tries to move second closer to first tensor
    if gamma_square_exists and (num1 != num2 and abs(num1 - num2) != 1): # won't trigger if we simplified gamma_squared
      # need to find spin_relation
      for relationA in term.relation_list:
        if relationA.index_type == "spinorial" and relationA.relation_type == "contraction":
          if (relationA.tensors[0] is gammaL and relationA.tensors[1] is gammaR) or (relationA.tensors[1] is gammaL and relationA.tensors[0] is gammaR):
            spin_relation = relationA
            break # found it

      apply_Clifford(expression, term, gammaL, gammaR, spin_relation)

# techincally trusty that C's, eta's reduced
def find_gamma_chains(term, mode = "only gammas", singletons = False, C = False):

  if mode == "only gammas":
    allowed_types = {"gamma"}
  if mode == "with gamma5s":
    allowed_types = {"gamma", "gamma5"}

  gamma_chains = []
  end = False
  while end == False:
    end = True # default
    for relation in term.relation_list:
      if relation.relation_type == "contraction" and relation.index_type == "spinorial":
        tensor1 = relation.tensors[0]
        tensor2 = relation.tensors[1]

        if tensor1.special_type in allowed_types and tensor2.special_type in allowed_types: # this is chain condition

          in_chain = False

          for i in range(len(gamma_chains)):
            if tensor1 in gamma_chains[i] and tensor2 in gamma_chains[i]:
              in_chain = True
              # no break statement so it goes faster through making duplicates
              # know that either tensor1 is on the end or tensor2 is on the end or both are in or niether are in
            elif tensor1 is gamma_chains[i][0]: # assuming no trace/gamma^2 anymore, so only one relation between gammas
              gamma_chains[i].insert(0, tensor2)   # want this to be adjacent the other one is
              in_chain = True
              end = False
              # no break statement so it goes faster through making duplicates
            elif tensor1 is gamma_chains[i][-1]: # can only be first or last in chain, middle ones already taken
              gamma_chains[i].append(tensor2)
              in_chain = True
              end = False
            elif tensor2 is gamma_chains[i][0]:
              gamma_chains[i].insert(0, tensor1)
              in_chain = True
              end = False
            elif tensor2 is gamma_chains[i][-1]:
              gamma_chains[i].append(tensor1)
              in_chain = True
              end = False

          if in_chain == False and tensor1 is not tensor2: # then we start a new chain
            gamma_chains.append([tensor1, tensor2])
            end = False
          elif in_chain == False and tensor1 is tensor2: # this is for self contractions for the trace
            gamma_chains.append([tensor1])
            end = False

  # made gamma_chains

  # now must delete repeat gamma_chain
  i = 0
  while i < len(gamma_chains):
    j = i + 1

    while j < len(gamma_chains):
      if gamma_chains[j][0] in gamma_chains[i]:
        del gamma_chains[j]
        j -= 1 # deindex because we removed one
      j += 1

    i += 1

  if singletons == True:
    gamma_chains_set = set()
    for gamma_chain in gamma_chains:
      for gamma in gamma_chain:
        gamma_chains_set.add(gamma)
    for tensor in term.tensor_list:
      if tensor.special_type in allowed_types and tensor not in gamma_chains_set:
        gamma_chains.append([tensor])

  if C == True:
    for tensor in term.tensor_list:
      if tensor.special_type == "C_matrix":
        gamma_chains.append([tensor])

  return gamma_chains

def TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list):
  # old_tensors will be a set of tensors in the term. got to remove these
  # old_relations will be a set of relations in the term. got to remove these
  # new_tensors_list is the list at index [0], this is the one to replace with in the term, the rest are for output terms
  # new_relations_list is similar as above
  # see loop for format of tensor_index_dit, these show which one to transfer to
  # will have to replace references to the old tensors with the new ones in other relations.

  output_terms = []

  # removing both old_tensors and old_relations early
  for tensor in old_tensors:
    term.tensor_list.remove(tensor)
  for relation in old_relations:
    term.relation_list.remove(relation)

  num_of_new_terms = len(new_tensors_list) - 1  # -1 because making use of original term too

  # creating new terms
  for i in range(num_of_new_terms):
    output_terms.append(term.copy())

  # reassigning relations for each term
  for i in range(len(term.relation_list)): # using the first term's relations to figure out how to change the rest
    relation = term.relation_list[i]
    for j in range(len(relation.tensors)):
      tensor = relation.tensors[j]
      index = relation.indices[j]
      if tensor in old_tensors: # we have a hit, can be a set for O(1) lookup
        # now find the new corresponding tensor
        for k in range(len(new_tensors_list)):
          if k == 0:
            relation.tensors[j] = tensor_index_dict[(tensor, index, k, relation.index_type)][0]
            relation.indices[j] = tensor_index_dict[(tensor, index, k, relation.index_type)][1]
          else:
            term_temp = output_terms[k-1]
            relation_temp = term_temp.relation_list[i]

            relation_temp.tensors[j] = tensor_index_dict[(tensor, index, k, relation.index_type)][0]
            relation_temp.indices[j] = tensor_index_dict[(tensor, index, k, relation.index_type)][1]


  # adding the rest of the tensors & relations we have to, order doesn't matter because will be reordered anyways and all indices explicit
  for i in range(len(new_tensors_list)):
    if i == 0:
      term.tensor_list += new_tensors_list[i]
      term.relation_list += new_relations_list[i]
      term.coefficient *= coeff_factor_list[i]
    else:
      output_terms[i-1].tensor_list += new_tensors_list[i] # output_terms[i-1] is ith term
      output_terms[i-1].relation_list += new_relations_list[i]
      output_terms[i-1].coefficient *= coeff_factor_list[i]

  return output_terms


def expand_all_sigmas(expression):
  i = 0
  while i < len(expression.term_list):
    term = expression.term_list[i]
    j = 0
    while j < len(term.tensor_list): # tensor_list will change, need while loop
      tensor = term.tensor_list[j]
      if tensor.special_type == "sigma":
        expand_sigma(expression, term, tensor)
        j -= 1 # deindexing since we removed the tensor
      j += 1
    i += 1

def expand_sigma(expression, term, tensor):
  # tensor is sigma
  sigma = tensor


  # create old_tensors set
  old_tensors = {sigma}

  # create old_relations set
  old_relations = set()

  # create new_tensors_list (has sublists per term)
  mu = tensor.index_structure[0][0]
  nu = tensor.index_structure[0][1]
  a = tensor.index_structure[1][0]
  b = tensor.index_structure[1][1]
  gammaL1 = Gamma_mu([[mu],[a,"up"]])
  gammaR1 = Gamma_mu([[nu],["down",b]])
  gammaL2 = Gamma_mu([[nu],[a,"up"]])
  gammaR2 = Gamma_mu([[mu],["down",b]])

  new_tensors_list = [[gammaL1,gammaR1],[gammaL2,gammaR2]]

  # create new_relations_list (has sublists per term)
  relation1 = Contraction([gammaL1,gammaR1],[1,0],"spinorial") # contraction has args (tensors, indices, index_type)
  relation2 = Contraction([gammaL2,gammaR2],[1,0],"spinorial")

  new_relations_list = [[relation1], [relation2]]

  # create tensor_index_dict
  tensor_index_dict = {}
  tensor_index_dict[(sigma, 0, 0, "spacetime")] = [gammaL1, 0] # tensor is sigma, tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(sigma, 1, 0, "spacetime")] = [gammaR1, 0]
  tensor_index_dict[(sigma, 0, 1, "spacetime")] = [gammaR2, 0]
  tensor_index_dict[(sigma, 1, 1, "spacetime")] = [gammaL2, 0]

  tensor_index_dict[(sigma, 0, 0, "spinorial")] = [gammaL1, 0]
  tensor_index_dict[(sigma, 1, 0, "spinorial")] = [gammaR1, 1]
  tensor_index_dict[(sigma, 0, 1, "spinorial")] = [gammaL2, 0]
  tensor_index_dict[(sigma, 1, 1, "spinorial")] = [gammaR2, 1]

  # create coeff_factor_list
  coeff_factor_list = [1j/2,-1j/2] # second one gets multiplied by -1


  new_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
  expression.term_list = expression.term_list + new_terms


def expand_all_gamma5(expression):

  i = 0
  while i < len(expression.term_list):
    term = expression.term_list[i]
    j = 0
    while j < len(term.tensor_list): # tensor_list will change, need while loop
      tensor = term.tensor_list[j]
      if tensor.special_type == "gamma5":
        trusty_expand_gamma5(expression, term, tensor)
        j -= 1 # deindexing since we removed the tensor
      j += 1
    i += 1

def TERMLEVELexpand_all_gamma5(expression, term):

  j = 0
  while j < len(term.tensor_list): # tensor_list will change, need while loop
    tensor = term.tensor_list[j]
    if tensor.special_type == "gamma5":
      trusty_expand_gamma5(expression, term, tensor)
      j -= 1 # deindexing since we removed the tensor
    j += 1

# VERY trusty that gamma trace does not exist
def trusty_expand_gamma5(expression, term, gamma5):

  gamma_chains = find_gamma_chains(term)

  adjacent_gammas = []

  for relation in term.relation_list:
    if relation.index_type == "spinorial" and relation.relation_type == "contraction":
      for i in range(len(relation.tensors)):
        tensor1 = relation.tensors[0]
        tensor2 = relation.tensors[1]
        if tensor2 is gamma5:
          adjacent_gammas.append(tensor1)

  # getting more gammas in gamma_to_delete
  gamma_chainL = []
  gamma_chainR = []
  old_tensors_list = [gamma5]

  # first finding gamma_chainL and gamma_chainR. neither are garaunteed to exist, but if gamma_chainR exists => gamma_chainL was found
  for gamma in adjacent_gammas:
    for gamma_chain in gamma_chains:
      if gamma in gamma_chain:
        if gamma_chainL == []:
          gamma_chainL = gamma_chain
        elif gamma_chain is not gamma_chainL:
          gamma_chainR = gamma_chain
          break
    if gamma_chainR != []:
      break # found both

  gamma_chainL_set = set(gamma_chainL) # for faster lookups later
  gamma_chainR_set = set(gamma_chainR)

  # can fit up to 3
  end = False
  while len(old_tensors_list) < 4 and end == False: # was len(gammas_to_delete) < 3
    end = True # default

    for gamma in old_tensors_list:
      if gamma is not gamma5:
        if gamma in gamma_chainL_set:
          j = gamma_chainL.index(gamma)
          for i in range(len(gamma_chainL)):
            if abs(i - j) == 1 and gamma_chainL[i] not in old_tensors_list:
              old_tensors_list.insert(0, gamma_chainL[i])
              end = False # found one, continue
              break
        elif gamma in gamma_chainR_set:
          j = gamma_chainR.index(gamma)
          for i in range(len(gamma_chainR)):
            if abs(i - j) == 1 and gamma_chainR[i] not in old_tensors_list:
              old_tensors_list.append(gamma_chainR[i])
              end = False
              break
        if end == False:
          break # go to the next while

  # create old_tensors set
  old_tensors = set(old_tensors_list)

  # create old_relations set
  old_relations = set()


  a_tensor = old_tensors_list[0]
  b_tensor = old_tensors_list[-1]

  if a_tensor is b_tensor: # this is the case for one gamma5
    a_index = 0
    b_index = 1
  else: # gamma5 plus other stuff
    hits = 0

    for relation in term.relation_list:
      if relation.index_type == "spinorial":
        if relation.relation_type == "contraction":
          for i in range(len(relation.tensors)):
            tensor1 = relation.tensors[0]
            tensor2 = relation.tensors[1]
            if tensor1 is a_tensor and tensor2 not in old_tensors: # this is assumign gamma trace is reduced
              a_index = relation.indices[0]
            elif tensor1 is b_tensor and tensor2 not in old_tensors:
              b_index = relation.indices[0]
        else:
          if relation.tensors[0] is a_tensor:
            a_index = relation.indices[0]
          elif relation.tensors[0] is b_tensor:
            b_index = relation.indices[0]

        if relation.tensors[0] is a_tensor and relation.tensors[1] is b_tensor: # then it is looping
          hits += 1
          if hits == 2: # then double contraction gamma and gamma5
            a_index = relation.tensors[0]
            b_index = relation.tensors[1]
            break

  a = a_tensor.index_structure[1][a_index]
  b = b_tensor.index_structure[1][b_index]

  # found a,b, look for mu, nu, etc

  mu_list = []

  for tensor in old_tensors_list:
    if tensor is not gamma5:
      mu_list.append(tensor.index_structure[0][0])

  while len(mu_list) < 4: # completing mu_list
    mu_list.append("up")

  Levi = Levi_Civita([mu_list,[]]) # convention is levi_civita with all indices up is positive

  new_relations = []
  new_gammas = []

  # creating the gammas and their contractions to Levi
  for i in range(5 - len(old_tensors_list)):
    if i == 0:
      first_index = a
    else:
      first_index = "down"
    if i == 4 - len(old_tensors_list):
      second_index = b
    else:
      second_index = "up"

    gamma = Gamma_mu([["down"],[first_index,second_index]])
    new_gammas.append(gamma)
    # spacetime relation
    Levi_index = i + (len(old_tensors_list) - 1)
    new_relations.append(Contraction([Levi,gamma],[Levi_index,0],"spacetime"))

  new_tensors_list = [[Levi] + new_gammas]

  # now making spinorial contrctations between them
  for i in range(len(new_gammas) - 1):
    gamma1 = new_gammas[i]
    gamma2 = new_gammas[i+1]
    new_relations.append(Contraction([gamma1,gamma2],[1,0],"spinorial"))

  new_relations_list = [new_relations]

  # creating tensor_index_dict
  tensor_index_dict = {}

  tensor_index_dict[(a_tensor, a_index, 0, "spinorial")] = [new_gammas[0], 0] # tensor is gamma5, tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(b_tensor, b_index, 0, "spinorial")] = [new_gammas[-1], 1]

  i = 0
  for tensor in old_tensors_list:
    if tensor is not gamma5:
      tensor_index_dict[(tensor, 0, 0, "spacetime")] = [Levi, i]
      i += 1

  # create coeff_factor_list
  coeff_factor = 1j # N.B. the 1j present in the gamma5 formula
  all_term_coeff_factor = 1

  # fixing spinorial contraction convention
  for i in range(len(old_tensors_list) - 1):
    tensor1 = old_tensors_list[i]
    tensor2 = old_tensors_list[i+1]
    next = False
    for relation in term.relation_list:
      for j in range(len(relation.tensors)):
        if tensor1 is relation.tensors[j] and tensor2 is relation.tensors[1 - j]:
          if tensor1.index_structure[1][relation.tensors[j]] == "down":
            all_term_coeff_factor *= -1
          next = True
          break
      if next == True:
        break

  all_term_coeff_factor *= (-1)**len(gamma_chainR) # since we moved these to the left of gamma5 using gamma5 gamma mu = - gamma mu gamma5

  coeff_factor *= 1/24
  for i in range(len(old_tensors_list) - 1):
    coeff_factor *= (4-i)

  if len(old_tensors_list) > 1:
    coeff_factor *= -1 # per our identites

  coeff_factor_list = [coeff_factor*all_term_coeff_factor]

  # making the extra terms
  if len(old_tensors_list) == 3:
    eta = Eta([[mu_list[0]],[mu_list[1]]])
    new_gamma5 = Gamma5([a,b])

    new_tensors_list.append([eta,gamma5])
    new_relations_list.append([])

    i = 0
    for tensor in old_tensors_list:
      if tensor is not gamma5:
        tensor_index_dict[(tensor, 0, 1, "spacetime")] = [eta, i]
        i += 1
    tensor_index_dict[(a_tensor, a_index, 1, "spinorial")] = [gamma5, 0]
    tensor_index_dict[(b_tensor, b_index, 1, "spinorial")] = [gamma5, 1]

    coeff_factor_list.append(all_term_coeff_factor)

  if len(old_tensors_list) == 4:
    # extra term 1
    eta = Eta([[mu_list[0]],[mu_list[1]]])
    new_gamma = Gamma_mu([[mu_list[2]],[a,"up"]])
    new_gamma5 = Gamma5(["down",b])
    i = 0
    for tensor in old_tensors_list:
      if tensor is not gamma5:
        if i < 2:
          tensor_index_dict[(tensor, 0, 1, "spacetime")] = [eta, i]
        else:
          tensor_index_dict[(tensor, 0, 1, "spacetime")] = [new_gamma, 0]
        i += 1
    tensor_index_dict[(a_tensor, a_index, 1, "spinorial")] = [new_gamma, 0]
    tensor_index_dict[(b_tensor, b_index, 1, "spinorial")] = [new_gamma5, 1]
    new_tensors_list.append([eta,new_gamma,new_gamma5])
    new_relations_list.append([Contraction([new_gamma, new_gamma5],[1,0],"spinorial")])
    coeff_factor_list.append(all_term_coeff_factor)
    # extra term 2
    eta = Eta([[mu_list[1]],[mu_list[2]]])
    new_gamma = Gamma_mu([[mu_list[0]],[a,"up"]])
    new_gamma5 = Gamma5(["down",b])
    i = 0
    for tensor in old_tensors_list:
      if tensor is not gamma5:
        if i == 1 or i == 2:
          tensor_index_dict[(tensor, 0, 2, "spacetime")] = [eta, i-1] # i/2 will be 0 or 1
        else:
          tensor_index_dict[(tensor, 0, 2, "spacetime")] = [new_gamma, 0]
        i += 1
    tensor_index_dict[(a_tensor, a_index, 2, "spinorial")] = [new_gamma, 0]
    tensor_index_dict[(b_tensor, b_index, 2, "spinorial")] = [new_gamma5, 1]
    new_tensors_list.append([eta,new_gamma,new_gamma5])
    new_relations_list.append([Contraction([new_gamma, new_gamma5],[1,0],"spinorial")])
    coeff_factor_list.append(all_term_coeff_factor)
    # extra term 3
    eta = Eta([[mu_list[0]],[mu_list[2]]])
    new_gamma = Gamma_mu([[mu_list[1]],[a,"up"]])
    new_gamma5 = Gamma5(["down",b])
    i = 0
    for tensor in old_tensors_list:
      if tensor is not gamma5:
        if i == 0 or i == 2:
          tensor_index_dict[(tensor, 0, 3, "spacetime")] = [eta, int(i/2)] # i/2 will be 0 or 1
        else:
          tensor_index_dict[(tensor, 0, 3, "spacetime")] = [new_gamma, 0]
        i += 1
    tensor_index_dict[(a_tensor, a_index, 3, "spinorial")] = [new_gamma, 0]
    tensor_index_dict[(b_tensor, b_index, 3, "spinorial")] = [new_gamma5, 1]
    new_tensors_list.append([eta,new_gamma,new_gamma5])
    new_relations_list.append([Contraction([new_gamma, new_gamma5],[1,0],"spinorial")])
    coeff_factor_list.append(-1 * all_term_coeff_factor) # N.B. -1

  new_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

  expression.term_list = expression.term_list + new_terms


def apply_Clifford(expression, term, tensor1, tensor2, relation):
  # tensor1 is first gamma
  # tensor2 is second gamma
  # relation is their spinorial contraction

  # create old_tensors set
  old_tensors = {tensor1,tensor2}

  # create old_relations set
  old_relations = {relation}

  # create new_tensors_list (has sublists per term)
  # identify which contracted_indices
  for i in range(len(relation.tensors)):
    if relation.tensors[i] is tensor1:
      contracted_index1 = relation.indices[i]
      contracted_index2 = relation.indices[1 - i] # flips 0 and 1
      break # found it

  if tensor1.special_type == "gamma":
    mu = tensor1.index_structure[0][0]
  if tensor2.special_type == "gamma":
    nu = tensor2.index_structure[0][0]
  a = tensor1.index_structure[1][1 - contracted_index1] # flips 0 and 1
  c1 = tensor1.index_structure[1][contracted_index1]
  c2 = tensor2.index_structure[1][contracted_index2]
  b = tensor2.index_structure[1][1 - contracted_index2] # flips 0 and 1

  if tensor2.special_type == "gamma":
    gamma_likeL = Gamma_mu([[nu],[a,c1]]) # means we don't have to account for sign switching height of c1, c2
  else:
    gamma_likeL = Gamma5([[],[a,c1]])
  if tensor1.special_type == "gamma":
    gamma_likeR = Gamma_mu([[mu],[c2,b]])
  else:
    gamma_likeR = Gamma5([[],[c2,b]])

  if c1 == "up":
    change_height_factor = 1 # only does it for the C
  else:
    change_height_factor = -1

  # create new_relations_list (has sublists per term)
  relation1 = Contraction([gamma_likeL,gamma_likeR],[1,0],"spinorial") # contraction has args (tensors, indices, index_type)

  new_relations_list = [[relation1]] # no relations for the second term

  # create tensor_index_dict
  tensor_index_dict = {}
  if tensor1.special_type == "gamma":
    tensor_index_dict[(tensor1, 0, 0, "spacetime")] = [gamma_likeR, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  if tensor2.special_type == "gamma":
    tensor_index_dict[(tensor2, 0, 0, "spacetime")] = [gamma_likeL, 0]
  tensor_index_dict[(tensor2, 1 - contracted_index2, 0, "spinorial")] = [gamma_likeR, 1]
  tensor_index_dict[(tensor1, 1 - contracted_index1, 0, "spinorial")] = [gamma_likeL, 0] # anything pointing to a still points to a
  # if they are gamma5's then we don't have to reassign indices

  if tensor1.special_type == "gamma5" and tensor2.special_type == "gamma5":
    coeff_factor_list = [1]
  else:
    coeff_factor_list = [-1] # this will be true for gamma5 gamma mu and gamma mu gamma nu

  if tensor1.special_type == "gamma5":
    if contracted_index1 == 0:         # getting gamma5 into matrix mult form
      coeff_factor_list[0] *= -1
  if tensor2.special_type == "gamma5":
    if contracted_index2 == 1:
      coeff_factor_list[0] *= -1

  new_tensors_list = [[gamma_likeL,gamma_likeR]]
  if tensor1.special_type == "gamma" and tensor2.special_type == "gamma":
    eta = Eta([[mu,nu],[]])
    C = C_matrix([[],[a,b]]) # no considerations necessary for this as to whether it will make the term negative
    new_tensors_list.append([eta,C])
    new_relations_list.append([])
    tensor_index_dict[(tensor1, 1 - contracted_index1, 1, "spinorial")] = [C, 0]
    tensor_index_dict[(tensor2, 1 - contracted_index2, 1, "spinorial")] = [C, 1]
    tensor_index_dict[(tensor1, 0, 1, "spacetime")] = [eta, 0]
    tensor_index_dict[(tensor2, 0, 1, "spacetime")] = [eta, 1]
    coeff_factor_list.append(2 * change_height_factor) # N.B. change_height_factor

  new_terms = TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
  expression.term_list = expression.term_list + new_terms


# trusts there are no etas left over
# trusts that levis have been dealt with
def trusty_auto_kill_terms(expression):
  i = 0
  while i < len(expression.term_list): # we will be removing terms
    term = expression.term_list[i]
    if TERMLEVELtrusty_auto_kill_terms(term):
      expression.term_list.remove(term)
      i -= 1 # deindexing because we just removed this term
    i += 1


# defining a function to check symmetry_structure
def symmetry(tensor, index1, index2, index_type):
  if index1 == index2:
    return None
  i = 0 if index_type == "spacetime" else 1
  if tensor.symmetry_structure[i] is not None and tensor.symmetry_structure[i] != []:
    for temp in tensor.symmetry_structure[i]:
      if (temp[0] == index1 and temp[1] == index2) or (temp[1] == index1 and temp[0] == index2):
        return temp[2]
  return None

def TERMLEVELtrusty_auto_kill_terms(term):

  TERMLEVELtrusty_reduce_eta_and_C(term)

  zeroed = False # if return True, then delete the term

  # even a single conraction to a levi with itself will kill it

  for i in range(len(term.relation_list)):
    relation1 = term.relation_list[i]
    if relation1.relation_type == "contraction":

      # Levi Civita self contraction case
      if relation1.tensors[0].special_type == "Levi Civita" and relation1.tensors[1] is relation1.tensors[0]: # will check all relations, don't worry about tensor2
        return True # have to do it here, in case there are no other contractions
        # using trusty that Levis reduced already
      # field self contraction case
      if relation1.tensors[0].special_type[:5] == "field" and relation1.tensors[1] is relation1.tensors[0]:
        # this has to be if there is a antisym spacetime field contracted to itself
        field = relation1.tensors[0]
        index1 = relation1.indices[0]
        index2 = relation2.indices[1]
        if relation1.relation_type == "spinorial" and "sym" == symmetry(field, index1, index2, "spinorial"): # special function defined at start

          # this has to be if there is a sym spin field (because index height flipping)
          return True

        if relation1.relation_type == "spacetime" and "antisym" == symmetry(field, index1, index2, "spacetime"):

          # this has to be if there is a antisym space field
          return True

      for j in range(i+1, len(term.relation_list)): # pairwise checking, no repeats or self checks
        relation2 = term.relation_list[j]
        if relation2.relation_type == "contraction" and relation1.index_type == relation2.index_type: # have to be same index_type for double contraction

          for k in range(2): # all mixes and matches
            tensor1A = relation1.tensors[k]
            tensor1B = relation1.tensors[1-k]
            for l in range(2):
              tensor2A = relation2.tensors[l]
              tensor2B = relation2.tensors[1-l]

              # this is a double contraction
              if (tensor1A is tensor2A and tensor1B is tensor2B):
                # then it is a double contraction, our double loop before checks all combos

                # only field zeroing left, since no contracted etas, C's left
                if tensor1A.special_type[:5] == "field":

                  field = tensor1A
                  index1 = relation1.indices[k]
                  index2 = relation2.indices[l]

                  if relation1.relation_type == "spinorial" and tensor1B.special_type == "gamma" and "sym" == symmetry(field, index1, index2, "spinorial"):

                    # has to be antisym field and contracts with gammamu ab
                    return True

                  if relation1.relation_type == "spacetime" and tensor1B.special_type == "Levi Civita" and ("sym" == symmetry(field, index1, index2, "spacetime") or "antisym" == symmetry(field, index1, index2, "spacetime")):
                    # either sym or antisym structure works, kills with Levi

                    # has to be sym or antisym field and contracts with Levi Civita
                    return True

              elif tensor1A.special_type == "partial" and tensor2A.special_type == "partial" and tensor1B is tensor2B:
                index1 = relation1.indices[k]
                index2 = relation2.indices[l]
                if tensor1B.special_type == "Levi Civita" or "antisym" == symmetry(tensor1B, index1, index2, "spacetime"):
                  return True

  return zeroed


# handler for converting connected gammas double contracted with partials into etas

# going to put this inside of the termwise steps
def reduce_two_gamma_two_partial(expression, term):

  end = False

  while end == False:
    end = True

    TERMLEVELtrusty_reduce_eta_and_C(term)

    gamma_chains = find_gamma_chains(term)

    hits = 0 # this is needed to instantiate it in case we have no gamma chains (for the if hits < 2)

    for gamma_chain_temp in gamma_chains:

      hits = 0
      gammas_to_reduce = []
      partials_to_reduce = []

      for relation in term.relation_list:
        if relation.index_type == "spacetime" and relation.relation_type == "contraction":
          for i in range(len(relation.tensors)):
            tensor1 = relation.tensors[i]
            tensor2 = relation.tensors[1-i]
            if tensor1 in gamma_chain_temp and tensor2.special_type == "partial":
              hits += 1
              gammas_to_reduce.append(tensor1)
              partials_to_reduce.append(tensor2)
              gamma_chain = gamma_chain_temp
              break
          if hits == 2:
            end = False
            break

    if hits < 2:
      break

    # now we get the gammas in gamma_to_reduce to be adjacent
    temp1 = gamma_chain.index(gammas_to_reduce[0])
    temp2 = gamma_chain.index(gammas_to_reduce[1])
    index_pos1 = min(temp1, temp2)
    index_pos2 = max(temp1, temp2)

    for swap_index in range(index_pos2, index_pos1 + 1, -1):
        gammaL = gamma_chain[swap_index-1]
        gammaR = gamma_chain[swap_index]
        for relation in term.relation_list:
          if relation.relation_type == "contraction" and relation.index_type == "spinorial":
            if (relation.tensors[0] is gammaL and relation.tensors[1] is gammaR) or (relation.tensors[1] is gammaL and relation.tensors[0] is gammaR): # trusting no space contraction between these as well
              relationA = relation
              break # found it
        apply_Clifford(expression, term, gammaL, gammaR, relationA)
        gamma_chain[swap_index] = term.tensor_list[-1] # gamma_chain[swap_index-1]
        gamma_chain[swap_index-1] = term.tensor_list[-2] # the second to last appended tensor will be here in the chain

    # now we can apply the idetity. partials will be the same but now source gammas from gamma chain
    gamma1 = gamma_chain[index_pos1]
    gamma2 = gamma_chain[index_pos1 + 1]
    partial1 = partials_to_reduce[0]
    partial2 = partials_to_reduce[1]

    #old_tensors = {gamma1, gamma2}
    old_tensors = {gamma1, gamma2, partial1, partial2}

    # gather all the relations between each other and with the partials
    old_relations = set()
    for relation in term.relation_list:
      if relation.relation_type == "contraction":
        if relation.tensors[0] in old_tensors and relation.tensors[1] in old_tensors:
          old_relations.add(relation)
          if relation.tensors[0] is gamma1 and relation.tensors[1] is gamma2:
            index1 = 1 - relation.indices[0]
            index2 = 1 - relation.indices[1]
          elif relation.tensors[1] is gamma1 and relation.tensors[0] is gamma2:
            index1 = 1 - relation.indices[1]
            index2 = 1 - relation.indices[0]

    # creating new tensors
    a = gamma1.index_structure[1][index1] # flips 0 and 1
    c1 = gamma1.index_structure[1][1 - index1]
    b = gamma2.index_structure[1][index2] # flips 0 and 1

    C = C_matrix([[],[a,b]]) # no considerations necessary for this as to whether it will make the term negative
    new_partial1 = Partial([["up"],[]])
    new_partial2 = Partial([["down"],[]])

    if c1 == "up":
      change_height_factor = 1 # only does it for the C
    else:
      change_height_factor = -1

    new_tensors_list = [[C,new_partial1,new_partial2]]

    # create new_relations_list (has sublists per term)
    relationA = Contraction([new_partial1,new_partial2],[0,0],"spacetime") # contraction has args (tensors, indices, index_type)

    new_relations_list = [[relationA]] # no relations for the second term


    # creating tensor_index_dict
    tensor_index_dict = {}

    tensor_index_dict[(gamma1, index1, 0, "spinorial")] = [C, 0]
    tensor_index_dict[(gamma2, index2, 0, "spinorial")] = [C, 1]

    # creating coeff_factor_list
    coeff_factor_list = [change_height_factor]


    TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)


# trusty because assumes gamma^2 has been dealt with, although don't think it matters the order
# and trusts no C's left
def trusty_gamma_trace(expression):
  i = 0
  while i < len(expression.term_list): # we will be removing terms
    term = expression.term_list[i]

    output_terms, zeroed = TERMLEVELtrusty_gamma_trace(term)

    if zeroed:
      expression.term_list.remove(term)
      i -= 1 # deindexing because we removed a term
    else:
      expression.term_list = expression.term_list + output_terms

    i += 1

def TERMLEVELrecursive_trace_formula(term, gamma_chain, output_terms = None):
  if output_terms is None:
    output_terms = []

  # building up args to use TERMLEVELapply_identity
  old_tensors = set(gamma_chain)
  old_relations = set()

  # let's identify old_relations for the gamma chain
  for relation in term.relation_list:
    if relation.tensors[0] in old_tensors and relation.index_type == "spinorial" and relation.relation_type == "contraction":
      old_relations.add(relation)  # know it's a trace so all spinorial gamma chain contractions will be there

  # end condition
  if len(gamma_chain) == 2:
    mu = gamma_chain[0].index_structure[0][0]
    nu = gamma_chain[1].index_structure[0][0]
    eta = Eta([[mu, nu],[]])

    new_tensors_list = [[eta]]

    new_relations_list = [[]] # none to append

    tensor_index_dict = {}
    tensor_index_dict[(gamma_chain[0], 0, 0, "spacetime")] = [eta, 0] # tuple is (tensor, index, term_num, relation.relation_type)
    tensor_index_dict[(gamma_chain[1], 0, 0, "spacetime")] = [eta, 1]

    if gamma_chain[0].index_structure[1][0] == gamma_chain[0].index_structure[1][1]: # this is if both indices are up/down, so it's not matrix multiplication
      change_height_factor = -1
    else:
      change_height_factor = 1

    coeff_factor_list = [change_height_factor * 4] # this is convention

    output_terms += TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)
    return output_terms # returns here, no need to recursion

  # if we are at this part of the code, length is still >2

  new_tensors_list = []
  new_relations_list = []
  tensor_index_dict = {}
  coeff_factor_list = []

  # each k will be a new term
  for k in range(1, len(gamma_chain)): # skipping the first one per the formula. Recall 0 indexing in python
    new_tensors = []
    new_relations = []

    # change_height_factor if the one at k is not matrix mult form
    if gamma_chain[k].index_structure[1][0] == gamma_chain[k].index_structure[1][1]: # this is if both indices are up/down, so its not matrix multiplication
      change_height_factor = -1
    else:
      change_height_factor = 1

    for i in range(1, len(gamma_chain)): # skipping first one per formula
      if i != k:

        # testing for i adjacent to k, with wrapping
        if i == k-1: # left tensor of kth
          index_structure = [list(gamma_chain[i].index_structure[0]), list(gamma_chain[i].index_structure[1])] # deep copy
          index_structure[1][1] = "up"
          gamma = Gamma_mu(index_structure)
        elif k == 1 and i == len(gamma_chain)-1: # left tensor of kth, wrapping
          index_structure = [list(gamma_chain[i].index_structure[0]), list(gamma_chain[i].index_structure[1])]
          index_structure[1][1] = "up"
          gamma = Gamma_mu(index_structure)

        elif i == k+1: # right tensor of kth
          index_structure = [list(gamma_chain[i].index_structure[0]), list(gamma_chain[i].index_structure[1])]
          index_structure[1][0] = "down"
          gamma = Gamma_mu(index_structure)
        elif k == len(gamma_chain)-1 and i == 1: # right tensor of kth, wrapping
          index_structure = [list(gamma_chain[i].index_structure[0]), list(gamma_chain[i].index_structure[1])]
          index_structure[1][0] = "down"
          gamma = Gamma_mu(index_structure)

        else: # there was no adjacency
          gamma = gamma_chain[i].copy()

        if k == 1: # have to go to right side of gamma_chain
          if i == 2:
            gamma1 = gamma_chain[-1] # have to wrap around
            gamma2 = gamma_chain[i]
          else:
            gamma1 = gamma_chain[i-1]
            gamma2 = gamma_chain[i]
        elif k == len(gamma_chain)-1: # have to go to left side of gamma_chain
          if i == len(gamma_chain)-2:
            gamma1 = gamma_chain[i]
            gamma2 = gamma_chain[1]
          else:
            gamma1 = gamma_chain[i]
            gamma2 = gamma_chain[i+1]

        new_tensors.append(gamma)

    if k == 1:
      gamma_chainA = list(new_tensors) # no eta yet, this should be ordered too

    coeff_factor_list.append(change_height_factor * (-1)**(k+1)) # +1 since index at 0

    for i in range(len(new_tensors)): # will only have gammas at this point, eta not spawned yet, also will be ordered
      if i == 0:
        gamma1 = new_tensors[-1] # have to wrap around
        gamma2 = new_tensors[0]
      else:
        gamma1 = new_tensors[i-1]
        gamma2 = new_tensors[i]
      relation = Contraction([gamma1,gamma2],[1,0],"spinorial")
      new_relations.append(relation)

    # also spawn eta, add to new_tensors after tensor_index_dict
    mu = gamma_chain[0].index_structure[0][0]
    nu = gamma_chain[k].index_structure[0][0]
    eta = Eta([[mu, nu],[]])

    # make tensor_index_dict now, before we append eta to new_tensors
    temp = 0
    for i in range(1, len(gamma_chain)):
      if i != k:
        tensor_index_dict[(gamma_chain[i], 0, k-1, "spacetime")] = [new_tensors[temp], 0] # tuple is (tensor, index, term_num, relation.relation_type)
        temp += 1
    # eta tensor_index_dict
    tensor_index_dict[(gamma_chain[0], 0, k-1, "spacetime")] = [eta, 0] # tuple is (tensor, index, term_num, relation.relation_type)
    tensor_index_dict[(gamma_chain[k], 0, k-1, "spacetime")] = [eta, 1]

    # k is the k-1th term

    new_tensors.append(eta)

    new_tensors_list.append(new_tensors)
    new_relations_list.append(new_relations)

  output_terms += TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

  return TERMLEVELrecursive_trace_formula(term, gamma_chainA, output_terms)


def TERMLEVELtrusty_gamma_trace(term):

  output_terms = []

  zeroed = False # if an odd chain trace it will return true

  gamma_chains = find_gamma_chains(term) # this is gamma_chains, self contracted gammas will be in here

  trace_list = []

  # test if the the last element is contracted to the first
  for gamma_chain in gamma_chains:
    for relation in term.relation_list:
      if relation.index_type == "spinorial" and relation.relation_type == "contraction":
        if (relation.tensors[0] is gamma_chain[0] and relation.tensors[1] is gamma_chain[-1]) or (relation.tensors[1] is gamma_chain[0] and relation.tensors[0] is gamma_chain[-1]):
          if len(gamma_chain) != 2:
            # these are the gamma_chains that are traces:
            # avoiding 2 since this if statement will always be true as one hit is sufficient
            # this one also covers a single gamma self contracted
            trace_list.append(gamma_chain)
            break
          elif len(gamma_chain) == 2:
            # we check if there is a second spinorial contraction that is not the first one

            for relation2 in term.relation_list:
              if relation2.index_type == "spinorial" and relation2.relation_type == "contraction":
                if (relation2 is not relation) and (relation2.tensors[0] is gamma_chain[0] and relation2.tensors[1] is gamma_chain[-1]) or (relation2.tensors[1] is gamma_chain[0] and relation2.tensors[0] is gamma_chain[-1]):
                  break # found it
                    # now len == 2 case is dealt with
            break # then there is no second contraction

  # trace_list made

  for gamma_chain in trace_list:
    if len(gamma_chain) % 2 == 1: # then there are an odd number of gammas, so it is 0
      return [], True # zeroed = True, no output_terms
    else:
      # now we apply the recursion relation
      # use TERMLEVELapply_identity
      output_terms += TERMLEVELrecursive_trace_formula(term, gamma_chain)

  return output_terms, zeroed


def reduce_Levi_Civita_products(expression):
  i = 0
  while i < len(expression.term_list): # we will be removing terms
    term = expression.term_list[i]
    output_terms = TERMLEVELreduce_Levi_Civita_products(term)
    expression.term_list = expression.term_list + output_terms
    i += 1


def TERMLEVELreduce_Levi_Civita_products(term):
  output_terms = []

  permutations = [
    [1, 2, 3, 4], [1, 2, 4, 3],
    [1, 3, 4, 2], [1, 3, 2, 4],
    [1, 4, 2, 3], [1, 4, 3, 2],
    [2, 1, 4, 3], [2, 1, 3, 4],
    [2, 3, 1, 4], [2, 3, 4, 1],
    [2, 4, 3, 1], [2, 4, 1, 3],
    [3, 1, 2, 4], [3, 1, 4, 2],
    [3, 2, 4, 1], [3, 2, 1, 4],
    [3, 4, 1, 2], [3, 4, 2, 1],
    [4, 1, 3, 2], [4, 1, 2, 3],
    [4, 2, 1, 3], [4, 2, 3, 1],
    [4, 3, 2, 1], [4, 3, 1, 2]
]

  # odd indexed permutations are odd, even are even (indexing at 0)
  # used for generating levi product terms

  end = False
  while end == False:
    end = True # default

    num_of_Levis = 0

    for tensor in term.tensor_list:
      if tensor.special_type == "Levi Civita":
        num_of_Levis += 1
        if num_of_Levis == 1:
          Levi1 = tensor
        if num_of_Levis == 2:
          Levi2 = tensor
          break
    if num_of_Levis < 2:
      break # no products of Levi's on this term!

    # if we are this part of the code, did not return
    end = False

    overall_factor = -1

    old_tensors = {Levi1, Levi2}
    old_relations = set()

    new_tensors_list = []
    new_relations_list = [] # no new relations

    muA = Levi1.index_structure[0][0]
    muB = Levi1.index_structure[0][1]
    muC = Levi1.index_structure[0][2]
    muD = Levi1.index_structure[0][3]
    nuA = Levi2.index_structure[0][0]
    nuB = Levi2.index_structure[0][1]
    nuC = Levi2.index_structure[0][2]
    nuD = Levi2.index_structure[0][3]

    nu_list = [nuA,nuB,nuC,nuD]

    tensor_index_dict = {}

    coeff_factor_list = []

    for term_num in range(len(permutations)):
      perm = permutations[term_num]
      i = perm[0]-1
      j = perm[1]-1
      k = perm[2]-1
      l = perm[3]-1

      nu1 = nu_list[i]
      nu2 = nu_list[j]
      nu3 = nu_list[k]
      nu4 = nu_list[l]

      eta1 = Eta([[muA, nu1],[]])
      eta2 = Eta([[muB, nu2],[]])
      eta3 = Eta([[muC, nu3],[]])
      eta4 = Eta([[muD, nu4],[]])

      new_tensors = [eta1, eta2, eta3, eta4]
      new_tensors_list.append(new_tensors)
      new_relations_list.append([]) # no new relations

      tensor_index_dict[(Levi1, 0, term_num, "spacetime")] = [eta1, 0]
      tensor_index_dict[(Levi1, 1, term_num, "spacetime")] = [eta2, 0]
      tensor_index_dict[(Levi1, 2, term_num, "spacetime")] = [eta3, 0]
      tensor_index_dict[(Levi1, 3, term_num, "spacetime")] = [eta4, 0]

      tensor_index_dict[(Levi2, i, term_num, "spacetime")] = [eta1, 1]
      tensor_index_dict[(Levi2, j, term_num, "spacetime")] = [eta2, 1]
      tensor_index_dict[(Levi2, k, term_num, "spacetime")] = [eta3, 1]
      tensor_index_dict[(Levi2, l, term_num, "spacetime")] = [eta4, 1]

      perm_sign = (-1)**term_num

      coeff_factor_list.append(overall_factor * perm_sign) # notice the overall_factor

    output_terms += TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

  return output_terms

# trusty because it assumes there is only one Levi in the term at this point
def trusty_Levi_Civita_gamma_converter(expression, term):
  # always go smaller

  Levi = None

  for tensor in term.tensor_list:
    if tensor.special_type == "Levi Civita":
      Levi = tensor
      break # found it

  if Levi is None:
    return

  # if our Levi is fully contracted to a chain we are good, don't do anything
  Levi_contraction_list = []
  for relation in term.relation_list:
    if relation.relation_type == "contraction":
      if relation.tensors[0] is Levi and relation.tensors[1].special_type == "gamma":
        Levi_contraction_list.append([relation.tensors[1], relation, relation.indices[0]]) # saving the relation for old_relations later, saving index for nu to avoid
      if relation.tensors[1] is Levi and relation.tensors[0].special_type == "gamma":
        Levi_contraction_list.append([relation.tensors[0], relation, relation.indices[1]])

  if len(Levi_contraction_list) == 0:
    return # no contractions with gammas to use the identity on!

  gamma_chains = find_gamma_chains(term)
  if len(gamma_chains) == 0:
    return
  for gamma_chain in gamma_chains:
    hits = 0

    for Levi_contraction in Levi_contraction_list:
      if Levi_contraction[0] in gamma_chain: # [0] is the tensor, [1] is the relation
        hits += 1


    if hits == 4:
      if len(gamma_chain) == 4:
        return # we are done, Levi is fully contracted to a chain, and we can't get it off 4

    elif hits < 4:
      return # we are done, Levi is fully contracted to a chain

  # if we reach this point, did not return

  # get the gamma5 gammas adjacent

  gamma_num_list = []

  for relation in term.relation_list:
    if relation.relation_type == "contraction" and relation.index_type == "spacetime":
      for i in range(len(relation.tensors)):
        tensor1 = relation.tensors[i]
        tensor2 = relation.tensors[1-i]
        if tensor1 in gamma_chain and tensor2 is Levi:
          gamma_num_list.append(gamma_chain.index(tensor1))
          # could add break statement but it is not necessary since trusty

  # now we get the gammas adjacent to each other
  gamma_num_list.sort()

  # since we sort it, going sequentially through the gamma in order does not affect the gamma_num of the later gammas

  for i in range(1, len(gamma_num_list)): # skipping 0th entry since we are getting the other gammas adjacent to this one, leave 0 alone
    gamma_num = gamma_num_list[i]
    if gamma_num - i != gamma_num_list[0]: # thus it is NOT already in the correct place
      for swap_index in range(gamma_num, gamma_num_list[0] + i, -1):
        gammaL = gamma_chain[swap_index-1]
        gammaR = gamma_chain[swap_index]

        for relation in term.relation_list:
          if relation.relation_type == "contraction" and relation.index_type == "spinorial":
            if (relation.tensors[0] is gammaL and relation.tensors[1] is gammaR) or (relation.tensors[1] is gammaL and relation.tensors[0] is gammaR): # trusting no space contraction between these as well
              relationA = relation
              break # found it

        apply_Clifford(expression, term, gammaL, gammaR, relationA)

        gamma_chain[swap_index] = term.tensor_list[-1] # gamma_chain[swap_index-1]
        gamma_chain[swap_index-1] = term.tensor_list[-2] # the second to last appended tensor will be here in the chain

  # now pre processing is complete, the gammas are all adjacent, we can prepare for apply_identity

  if gamma_num_list[0] == 0: # then the 5th gamma with be on the right, use gamma mu gamma5 = - gamma5 gamma mu to switch
    coeff_factor = -4
    gamma_to_destroy = gamma_chain[gamma_num_list[0] + 4]
  else: # then we can use the left gamma no problem
    gamma_to_destroy = gamma_chain[gamma_num_list[0] - 1]
    coeff_factor = 4


  old_tensors = {gamma_chain[gamma_num_list[0]], gamma_chain[gamma_num_list[0] + 1],
                 gamma_chain[gamma_num_list[0] + 2], gamma_chain[gamma_num_list[0] + 3],
                 gamma_to_destroy, Levi}

  reorder_list = [] # to get reordering Levi sign factor

  old_relations = set() # have to destroy all the connections between the gammas and with the Levi
  for relation in term.relation_list:
    # first the intermediate gamma contractions
    if relation.relation_type == "contraction" and relation.index_type == "spinorial":

      if relation.tensors[0] in old_tensors and relation.tensors[1] in old_tensors:
        old_relations.add(relation) # no gamma traces exist at this point in the code

        # have to consider spinorial contraction convention here!

        if gamma_chain.index(relation.tensors[0]) < gamma_chain.index(relation.tensors[1]): # using ordering of gamma_chain to determine spinorial contraction convention
          gamma1 = relation.tensors[0]
          index1 = relation.indices[0]
        else:
          gamma1 = relation.tensors[1]
          index1 = relation.indices[1]

        if gamma1.index_structure[1][index1] == "up":
          coeff_factor *= -1 # cost of flipping this contraction's two heights

    # now the gamma contractions to the Levi
    if relation.relation_type == "contraction" and relation.index_type == "spacetime":
      for i in range(len(relation.tensors)):
        tensor1 = relation.tensors[i]
        tensor2 = relation.tensors[1-i]
        if tensor1 in old_tensors and tensor2 is Levi: # covers Levi self contractions too accidentally, although there will be none
          old_relations.add(relation) # trusty that there is one Levi at this point

          # also need to sign factor form reordering Levi
          reorder_list.append([relation.indices[1-i], tensor1]) # [Levi's index, gamma]

  # now fixing sign from Levi reordering pre application of the identity (turning into gamma5)
  reorder_list.sort(key=lambda x: x[0]) # sorting by indices
  for num in range(4):
    from_num = gamma_chain.index(reorder_list[i][1]) - gamma_num_list[0]
    if num != from_num:
      temp1 = reorder_list[num]
      temp2 = reorder_list[from_num]
      reorder_list[num] = temp2
      reorder_list[from_num] = temp1

      coeff_factor *= -1 # now it is fixed into a gamma5


  # new_tensors_list
  new_Levi = Levi.copy() # making it a copy
  new_Levi.index_structure[0][0] == gamma_to_destroy.index_structure[0][0]

  mu = "up" if new_Levi.index_structure[0][1] == "down" else "down"
  nu = "up" if new_Levi.index_structure[0][2] == "down" else "down"
  alpha = "up" if new_Levi.index_structure[0][3] == "down" else "down"
  # finding a and b (ends of gamma subchain)
  for relation in term.relation_list:
    if relation.relation_type == "contraction" and relation.index_type == "spinorial":
      for i in range(len(relation.tensors)):
        tensor1 = relation.tensors[i]
        tensor2 = relation.tensors[i-1]
        if gamma_num_list[0] == 0: # case where gamma_to_destory is on the right
          if tensor1 is gamma_to_destroy and tensor2 in old_tensors: # assuming no traces
            b_index = 1 - relation.indices[i] # 1- flips 0 and 1
            b_tensor = tensor1
          if tensor1 is gamma_chain[gamma_num_list[0]] and tensor2 in old_tensors:
            a_index = 1 - relation.indices[i]
            a_tensor = tensor1
        else: # case where gamma_to_destroy is on the left
          if tensor1 is gamma_to_destroy and tensor2 in old_tensors: # assuming no traces
            a_index = 1 - relation.indices[i] # 1- flips 0 and 1
            a_tensor = tensor1
          if tensor1 is gamma_chain[gamma_num_list[0] + 3] and tensor2 in old_tensors : # other extreme of the subchain
            b_index = 1 - relation.indices[i]
            b_tensor = tensor1
  # found a_index, b_index, a_tensor, b_tensor
  a = a_tensor.index_structure[1][a_index]
  b = b_tensor.index_structure[1][b_index]

  gamma1 = Gamma_mu([[mu],[a,"down"]])
  gamma2 = Gamma_mu([[nu],["up","down"]])
  gamma3 = Gamma_mu([[alpha],["up",b]])

  new_tensors_list = [[new_Levi, gamma1, gamma2, gamma3]]

  # new_relations_list
  relation1 = Contraction([gamma1,gamma2],[1,0],"spinorial")
  relation2 = Contraction([gamma2,gamma3],[1,0],"spinorial")

  relation3 = Contraction([new_Levi,gamma1], [1,0],"spacetime")
  relation4 = Contraction([new_Levi,gamma2], [2,0],"spacetime")
  relation5 = Contraction([new_Levi,gamma3], [3,0],"spacetime")

  new_relations_list = [[relation1,relation2,relation3,relation4,relation5]]

  # tensor_index_dict
  tensor_index_dict = {}

  tensor_index_dict[(a_tensor, a_index, 0, "spinorial")] = [gamma1, 0] # tuple is (tensor, index, term_num, relation.relation_type)
  tensor_index_dict[(b_tensor, b_index, 0, "spinorial")] = [gamma3, 1]
  tensor_index_dict[(gamma_to_destroy, 0, 0, "spacetime")] = [new_Levi, 0]


  # coeff_factor_list
  coeff_factor_list = [coeff_factor] # have been editing coeff_factor as necessary above



  # will be no output terms
  TERMLEVELapply_identity(term, old_tensors, old_relations, new_tensors_list, new_relations_list, tensor_index_dict, coeff_factor_list)

