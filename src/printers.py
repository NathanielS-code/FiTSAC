# file printers

from objects import *

from fractions import Fraction # importing this to make output prettier

# making a printer for the expression
def expression_to_LaTeX(expression):

  output_string = ""

  for i in range(len(expression.term_list)):
    term = expression.term_list[i]
    part_of_string = term_to_LaTeX(term)
    if i != 0:
      output_string += " + "
    output_string += part_of_string

  # outputs a LaTeX string
  return output_string


# printer for a term:
def term_to_LaTeX(term):

  space_list = [
    r"\alpha", r"\beta", r"\gamma", r"\delta", r"\epsilon", # these are raw strings, prevent unintended escapes in python from the backslashes
    r"\zeta", r"\eta", r"\theta", r"\iota", r"\kappa",
    r"\lambda", r"\mu", r"\nu", r"\xi", r"\omicron",
    r"\pi", r"\rho", r"\sigma", r"\tau", r"\upsilon",
    r"\phi", r"\chi", r"\psi", r"\omega",
    r"\Alpha", r"\Beta", r"\Gamma", r"\Delta", r"\Epsilon",
    r"\Zeta", r"\Eta", r"\Theta", r"\Iota", r"\Kappa",
    r"\Lambda", r"\Mu", r"\Nu", r"\Xi", r"\Omicron",
    r"\Pi", r"\Rho", r"\Sigma", r"\Tau", r"\Upsilon",
    r"\Phi", r"\Chi", r"\Psi", r"\Omega"
  ]

  spin_list = [
      "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z",
      "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"
  ]


  space_label_tracker = 0
  spin_label_tracker = 0


  special_type_dict = {
      "gamma" : r"\gamma",  # these are raw strings, prevent unintended escapes in python from the backslashes
      "gamma5" : r"\gamma^5",
      "sigma" : r"\sigma",
      "C_matrix" : "C",
      "metric" : r"\eta",
      "partial" : r"\partial",
      "Levi Civita" : r"\epsilon",
  }

  # make relation_dict now
  relation_dict = {}

  free_space = []
  free_spin = []

  for relation in term.relation_list: # just storing which free's are space and spin
    if relation.relation_type[:4] == "free":
      if relation.index_type == "spacetime":
        free_space.append(relation.relation_type)
      if relation.index_type == "spinorial":
        free_spin.append(relation.relation_type)

  free_space = sorted(free_space, key=lambda x: int(x[4:]))
  free_spin = sorted(free_spin, key=lambda x: int(x[4:]))

  # free relations
  for relation in term.relation_list:
    if relation.relation_type[:4] == "free":
      if relation.index_type == "spinorial":

        temp = free_spin.index(relation.relation_type)

        relation_dict[relation] = spin_list[temp] # assigns space_label
        spin_label_tracker += 1

      elif relation.index_type == "spacetime":

        temp = free_space.index(relation.relation_type)

        relation_dict[relation] = space_list[temp] # assigns space_label
        space_label_tracker += 1

  # now contractions
  for relation in term.relation_list:
    if relation.relation_type == "contraction":
      if relation.index_type == "spinorial":

        relation_dict[relation] = spin_list[spin_label_tracker] # assigns spin_label
        spin_label_tracker += 1

      elif relation.index_type == "spacetime":

        relation_dict[relation] = space_list[space_label_tracker] # assigns spin_label
        space_label_tracker += 1

  # made relation_dict now

  output_string = r""

  coeff = term.coefficient
  real_part = Fraction(coeff.real).limit_denominator()
  imag_part = Fraction(coeff.imag).limit_denominator()


  if real_part != 0 and imag_part != 0:
    output_string += "(("

    if real_part.denominator == 1:
      output_string += str(real_part.numerator) + ")"
    else:
      output_string += "\\frac{" + str(real_part.numerator) + r"}{" + str(real_part.denominator) + "})"

    output_string += " + i("

    if imag_part.denominator == 1:
      output_string += str(imag_part.numerator) + ")"
    else:
      output_string += "\\frac{" + str(imag_part.numerator) + r"}{" + str(imag_part.denominator) + "})"

    output_string += ")) "

  else:
    if imag_part == 0:
      num = real_part
    if real_part == 0:
      num = imag_part

    if num < 0:
      output_string += "("

    if num.denominator == 1:
      if num.numerator != 1:
        output_string += str(num.numerator)
    else:
      output_string += "\\frac{" + str(num.numerator) + r"}{" + str(num.denominator) + "}"

    if num < 0:
      output_string += ")"
    else:
      output_string += ""

    if real_part == 0:
      output_string += "i "
    else:
      output_string += " "

  for tensor in term.tensor_list:

    # writing the tensor
    if tensor.special_type in special_type_dict:
      if tensor.special_type[:5] != "field":
        output_string += special_type_dict[tensor.special_type]
    else:
      if tensor.special_type[:5] != "field":
        output_string += tensor.special_type
      else:
        output_string += tensor.special_type[5] + "^{(" + tensor.special_type[6:] + ")}"

    # now we write the indices
    space_string = ""
    spin_string = ""

    # doing space in this loop
    i = 0 # about to += 1
    while i < len(tensor.index_structure[0]):

      hit = False

      for relation in term.relation_list:
        for j in range(len(relation.tensors)):
          tensor2 = relation.tensors[j]
          if tensor2 is tensor and relation.indices[j] == i and relation.index_type == "spacetime": # checking that we have the right tensor and index
            if tensor.index_structure[0][i] == "up":
              space_string += "{}^" + relation_dict[relation] # add to space_string
              hit = True

            elif tensor.index_structure[0][i] == "down":
              space_string += "{}_" + relation_dict[relation] # add to space_string
              hit = True
            break # found it
      if hit == False:
        spin_string += "{}_{Error} "
      i += 1
    # now spin index loop
    i = 0 # about to += 1
    while i < len(tensor.index_structure[1]):

      hit = False

      for relation in term.relation_list:

        for j in range(len(relation.tensors)):
          tensor2 = relation.tensors[j]
          if tensor2 is tensor and relation.indices[j] == i and relation.index_type == "spinorial": # checking that we have the right tensor and index and index type
            if tensor.index_structure[1][i] == "up":
              spin_string += "{}^" + relation_dict[relation] # add to spin_string
              hit = True

            elif tensor.index_structure[1][i] == "down":
              spin_string += "{}_" + relation_dict[relation] # add to spin_string
              hit = True
            break # found it
      if hit == False:
        spin_string += "{}_{Error} "
      i += 1


    output_string += space_string + spin_string + " " # combine the strings


  # outputs a LaTeX string
  return output_string



