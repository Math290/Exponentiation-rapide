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
//Calculer le double d'un point en coordonnées Jacobiennes
Double_in_Jacobi:=function(P)
    x1:=P[1];
    y1:=P[2];
    z1:=P[3];
    A:=4*x1*y1^2;
    B:=3*x1^2+a*z1^4;
    x3:=-2*A+B^2;
    y3:=-8*y1^4+B*(A-x3);
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
    A:=x1*(z2^2);
    B:=Q[1]*(z1^2);
    C:=y1*(z2^3);
    D:=y2*(z1^3);
    E:=B-A;
    F:=D-C;
    x3:=-(E^3)-2*A*(E^2)+(F^2);
    y3:=-C*(E^3) + F*(A*(E^2)-x3);
    z3:=z1*z2*E;
    return [x3,y3,z3];
end function;
question4:=function(n,P)
    R:=[Fp!1,Fp!1,Fp!0];
    bin:=Intseq(n,2);
    Reverse(~bin);
    for k in bin do
        R:=Double_in_Jacobi(R);
        if k eq 1 then
            R:=Somme_Jacobi(R,P);
        end if;
    end for;
    return Coordonnees_affines(R,E);
end function;
B:=Random(E); n:=Random(r); n*B eq question4(n,B);
A:=G; time for i:=1 to 100 do n:=Random(r);A:=question4(n,A); end for;



