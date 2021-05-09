//------------------------Méthode BRIP-----------------------------
p :=6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151;
r := 6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449;
Fp :=GF(p) ;
b := Fp !1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984 ;
Gx := Fp !2661740802050217063228768716723360960729859168756973147706671368418802944996427808491545080627771902352094241225065558662157113545570916814161637315895999846 ;
Gy := Fp !3757180025770020463545507224491183603594455134769762486694567779615544477440556316691234405012945539562144444537289428522585666729196580810124344277578376784 ;
a := Fp !-3 ;
E :=EllipticCurve([a,b]) ;
G :=E ![Gx,Gy] ;
//Passage de coordonnées jacobiennes en coordonnées affines.
Coordonnees_affines:=function(P,E)
    x1:=P[1];
    y1:=P[2];
    z1:=P[3];
    if z1 eq Fp!0 then
        return E!0;
    end if;
    return(E![x1/(z1^2), y1/(z1^3)]);
end function;
//Passage de coordonnées jacobiennes en coordonnées affines.
Coordonnees_affines:=function(P,E)
    x1:=P[1];
    y1:=P[2];
    z1:=P[3];
    if z1 eq Fp!0 then
        return E!0;
    end if;
    z1carre:=z1*z1;
    Z1:=z1carre*z1;
    return(E![x1/(z1carre), y1/(Z1)]);
end function;
//Calculer le double d'un point en coordonnées Jacobiennes
Double_in_Jacobi:=function(P)
    x1:=P[1];
    y1:=P[2];
    z1:=P[3];
    y1carre:=y1*y1;
    z1carre:=z1*z1;
    x1carre:=x1*x1;
    L:=x1*y1carre;
    A:=4*L;
    T:=z1carre*z1carre;
    B:=3*(x1carre-T);
    Bcarre:=B*B;
    x3:=-2*A+Bcarre;
    y3:=-8*y1carre*y1carre+B*(A-x3);
    z3:=2*y1*z1;
    return [x3,y3,z3];
end function;
//Somme en coordonnées jacobiennes
Somme_Jacobi:=function(P,Q) 
    x1:=P[1];
    y1:=P[2];
    z1:=P[3];
    x2:=Q[1];
    y2:=Q[2];
    z2:=Q[3];
    //Point à l'infini
    if x1 eq 1 and y1 eq 1 and z1 eq 0 then
        return Q;
    end if;
    //Point à l'infini
    if x2 eq 1 and y2 eq 1 and z2 eq 0 then
        return P;
    end if;
    z1carre:=z1*z1;
    z2carre:=z2*z2;
    A:=x1*(z2carre);
    B:=Q[1]*(z1carre);
    C:=y1*(z2carre*z2);
    D:=y2*(z1carre*z1);
    E:=B-A;
    F:=D-C;
    Ecarre:=E*E;
    Fcarre:=F*F;
    x3:=-(Ecarre*E)-2*A*(Ecarre)+(Fcarre);
    y3:=-C*(Ecarre*E) + F*(A*(Ecarre)-x3);
    z3:=z1*z2*E;
    return [x3,y3,z3];
end function;
Precomputations:=function(P,n)
    D:=Double_in_Jacobi(P);
    temp:=P;
    L:=[temp];
    k:=2^(n-1);
    for i:=1 to k do
        temp:=Somme_Jacobi(temp,D);
        Append(~L,temp);
    end for;
    return L;
end function;
// La fonction modulo signée
Mods:=function(d,w)
    k:=2^w;
    l:=2^(w-1);
    if d mod k ge l then
        return (d mod k) - k;
    else
        return d mod k;
    end if;
end function;
form_non_adjacente:=function(d,w)
    bin:=[];
    i:=0;
    while d gt 0 do
        if (d mod 2) eq 1 then
            Append(~bin,Mods(d,w));
            d:=d-bin[i+1];
        else
            Append(~bin,0);
        end if;
        d:= d div 2;
        i:=i+1;
    end while;
    return bin;
end function;
question4:=function(k,P)
    R:=Random(E);
    T:=Somme_Jacobi(P,-R);
    Q:=R;
    L:=Reverse(Intseq(k,2));
    t:=#L;
    for i:=1 to t do
        Q:=Double_in_Jacobi(Q);
        if L[i] eq 1 then
            Q:=Somme_Jacobi(Q,T);
        else
            Q:=Somme_Jacobi(Q,-R);
        end if;
    end for;
    S:=Somme_Jacobi(Q,-R);
    return Coordonnees_affines(S,E);
end function;

B:=Random(E); n:=Random(r); n*B eq question4(n,B);
A:=G; time for i:=1 to 100 do n:=Random(r);A:=question4(n,A); end for;

