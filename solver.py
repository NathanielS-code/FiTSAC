from tools import *
from objects import *
from printers import *

# solver goes here

def check_closure(field_list, transformation_list, print_level=1):

  algebra_satisfied = True

  for i in range(len(field_list)):

    LHS, RHS = setup_SUSY_algebra(i, field_list, transformation_list)

    for term in RHS.term_list:
      term.coefficient *= -1 # we are subtracting LHS from RHS

    LHSminusRHS = Expression(RHS.term_list + LHS.term_list)

    # now we get in to solution path

    expression = LHSminusRHS

    # 1st algorithm:

    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1:
      print("initial:")
      print(new_str)
    prev_str = new_str

    expand_all_sigmas(expression)
    new_str = "0 = " + expression_to_LaTeX(expression)

    if print_level >= 1 and prev_str != new_str:
      print("expanded sigmas:")
      print(new_str)
    prev_str = new_str

    for term in expression.term_list:
      trusty_reduce_gamma_squared(expression, term)
    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1 and prev_str != new_str:
      print("reduced gamma^2:")
      print(new_str)
    prev_str = new_str

    reduce_Levi_Civita_products(expression)
    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1 and prev_str != new_str:
      print("reduced Levi Civita products:")
      print(new_str)
    prev_str = new_str

    combine_like_terms(expression, print_level, mode = "heavy")
    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1: # no if statement for new_str here
      print("combined like terms:")
      print(new_str)
    prev_str = new_str

    trusty_reduce_eta_and_C(expression)

    # 2nd algorithm:

    expand_all_gamma5(expression)
    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1 and prev_str != new_str:
      print("expanded all gamma5s:")
      print(new_str)
    prev_str = new_str

    reduce_Levi_Civita_products(expression)
    new_str = "0 = " + expression_to_LaTeX(expression)
    if print_level >= 1 and prev_str != new_str:
      print("reduced Levi Civita products:")
      print(new_str)
    prev_str = new_str

    combine_like_terms(expression, print_level, mode = "super heavy")

    # now we check leftover terms
    if len(LHSminusRHS.term_list) != 0:
      algebra_satisfied = False
      print("|| Algebra satisfied for field " + str(i + 1) + " up to:")
      print("0 = " + expression_to_LaTeX(LHSminusRHS))
    else:
      print("|| Algebra satisfied for field " + str(i + 1))

  if algebra_satisfied:
    print("## Algebra satisfied overall")
  else:
    print("## Algebra !NOT! immediately satisfied for these fields (up to gauge & on-shell terms)")

  return algebra_satisfied
