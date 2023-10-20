######################################### FODEProblem ############################################

ftoomany(u, p, t, x, y) = 2u
order = 1.5
u0 = 0.5
tspan = (0.0, 1.0)
@test_throws SciMLBase.TooManyArgumentsError FODEProblem(ftoomany, order, u0, tspan)

ftoofew(u, t) = 2u
@test_throws SciMLBase.TooFewArgumentsError FODEProblem(ftoofew, order, u0, tspan)

fmessedup(u, t) = 2u
fmessedup(u, p, t, x, y) = 2u
@test_throws SciMLBase.FunctionArgumentsError FODEProblem(fmessedup, order, u0, tspan)
