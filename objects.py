# object definitions go here

# index_structure: is a 2-tuple with two nested n-tuples, e.g. (('up','down','down'),()).
# The first n-tuple contains the up and down indices of the spacetime indices,
# the second n-tuple contains the up and down indices of the spinorial indices

# symmetry_structure: is a 2-tuple containing two lists that each contain n
# 3-tuples, e.g. ([(1,2,'sym')],[]) indicates that the second and third
# spacetime indices (indexing from 0) are symmetric.
# General form of 3-tuple: (index1, index2, symmetry type),
# where symmetry type can be either 'sym' or 'antisym'
# NOTE: if indices are different heights, this is to be handled manually in solver

from decimal import Decimal
from fractions import Fraction

class Tensor:
  def __init__(self, index_structure, symmetry_structure=([],[]), special_type=None):
    self.index_structure = index_structure
    self.symmetry_structure = symmetry_structure
    self.special_type = special_type

  def copy(self):
    space_part = list(self.index_structure[0])
    spin_part = list(self.index_structure[1])
    new_index_structure = [space_part, spin_part]
    return Tensor(new_index_structure, self.symmetry_structure, self.special_type)

  def equals(self, tensor):
    if tensor.index_structure == self.index_structure and tensor.symmetry_structure == self.symmetry_structure and tensor.special_type == self.special_type:
      return True
    else:
      return False

#special versions
class Eta(Tensor):  # Eta is used to represent deltas as well
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([(0,1,"sym")],[]), special_type="metric")
class C_matrix(Tensor):  # C_matrix is used to represent deltas as well
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([],[(0,1,"antisym")]), special_type="C_matrix")
class Levi_Civita(Tensor):
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([(0,1,"antisym"),(1,2,"antisym"),(2,3,"antisym"), (0,3,"antisym"),
                                                                            (0,2,"antisym"), (1,3,"antisym")],[]), special_type="Levi Civita")
class Gamma_mu(Tensor):
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([],[(0,1,"sym")]), special_type="gamma")
class Gamma5(Tensor):
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([],[(0,1,"antisym")]), special_type="gamma5")
class Sigma(Tensor):
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([],[(0,1,"sym")]), special_type="sigma")

# -------------------------------------------
# used for tagmaking, only a temporary object
class Gamma_Chain(Tensor):
  def __init__(self, index_structure, num_of_gammas, num_of_gamma5s):
    super().__init__(index_structure=index_structure, symmetry_structure="Throw Error", special_type="gamma_chain" + str(num_of_gammas) + ";" + str(num_of_gamma5s))
    # should not access symmetry_structure
# -------------------------------------------

class Partial(Tensor):
  def __init__(self, index_structure):
    super().__init__(index_structure=index_structure, symmetry_structure=([],[]), special_type="partial")

class Field(Tensor):
  def __init__(self, index_structure, symmetry_structure, field_type, ID):
    special_type = "field" + field_type + str(ID) # field_type = "F" for fermionic or "B" for bosonic
    super().__init__(index_structure=index_structure, symmetry_structure=symmetry_structure, special_type=special_type)

class Relation:
  def __init__(self, tensors, indices, index_type, relation_type="contraction"):
    self.tensors = tensors # a list, could be of length 1 for free index. can have the same tensor at multiple points in this list, if contracted with itself
    self.indices = indices # a list, each indices[i] is the index related to tensors[i]
    self.index_type = index_type # "spacetime" or "spinorial"
    self.relation_type = relation_type # this becomes ID if working with a Free Index

class Contraction(Relation):
  def __init__(self, tensors, indices, index_type):
    super().__init__(tensors=tensors, indices=indices, index_type=index_type, relation_type="contraction")

class Free_Index(Relation):
  def __init__(self, tensor, index, index_type, ID):
    super().__init__(tensors=[tensor], indices=[index], index_type=index_type, relation_type="free"+str(ID))

class Term:
  def __init__(self, coefficient, tensor_list, relation_list):
    self.coefficient = coefficient  # this is just a number, e.g. "5" or "-2i"
    self.tensor_list = tensor_list  # order matters for tensor_list
    self.relation_list = relation_list

  @property
  def coefficient(self):
      real, imag = self._coefficient
      return complex(float(real), float(imag))

  @coefficient.setter
  def coefficient(self, value):
      c = complex(value)  # handles '1j', '2+3j', etc.
      self._coefficient = (
          Fraction(c.real).limit_denominator(),
          Fraction(c.imag).limit_denominator()
      )

  def copy(self):
    new_tensor_list = []
    new_relation_list = []

    for tensor in self.tensor_list:
      new_tensor_list.append(tensor.copy())

    for relation in self.relation_list: # now we copy over the relations
      new_tensors = []
      for i in range(len(relation.tensors)):
        tensor = relation.tensors[i]

        num = -1
        for i, tensorA in enumerate(self.tensor_list):
            if tensorA is tensor:
                num = i
                break

        if num != -1: # tensor in self.tensor_list:
          num = self.tensor_list.index(tensor)
          new_tensors.append(new_tensor_list[num])
        else:
          new_tensors.append(tensor) # this is necessary for apply_identity

      new_relation_list.append(Relation(list(new_tensors), list(relation.indices), relation.index_type, relation.relation_type))

    return Term(self.coefficient, new_tensor_list, new_relation_list)


# Expression
class Expression:
  def __init__(self, term_list):
    self.term_list = term_list

  def copy(self):
    new_term_list = []
    for term in self.term_list:
      new_term_list.append(term.copy())
    return Expression(new_term_list)
