using FractionalDiffEq

prob = QiSystem([35, 8/3, 80, -1, 1], [0.98, 0.98, 0.98], 50, [0.1, 0.2, 0.3])

solve(prob, QiAlg())