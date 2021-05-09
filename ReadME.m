----------------------------------------------------Version 5------------------------------------------------------------------
Pour cette version j'ai tenté d'optimiser la version précédente (W-NAF). 
Dans la fonction de doublement on peut changer B = 3X_{1}^2+aZ_{1}^4 par B=3(X_{1}^2-Z_{1}^4) puisque a=-3
D'autre part on remplace tous les squares par une multiplication (exemple A^2=A*A) et on précalcule 2^w et 2^(w-1).
Ce qui m'a permis de passer de 1.22 secondes à 1.06 secondes.
