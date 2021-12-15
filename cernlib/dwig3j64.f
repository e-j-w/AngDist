*
* $Id: dwig3j64.F,v 1.1.1.1 1996/04/01 15:01:48 mclareni Exp $
*
* $Log: dwig3j64.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:48  mclareni
* Mathlib gen
*
*
      FUNCTION DWIG3J(A1,B1,C1,X1,Y1,Z1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCG,LJN,LRC
 
      DIMENSION U(0:202)
 
      PARAMETER (R1 = 1, HF = R1/2)
 
      DATA U(0),U(2),(U(2*N-1),N=1,101) /103*0/
      DATA U(  4),U(  6) /6.931471805599453D-01, 1.791759469228055D+00/
      DATA U(  8),U( 10) /3.178053830347946D+00, 4.787491742782046D+00/
      DATA U( 12),U( 14) /6.579251212010101D+00, 8.525161361065414D+00/
      DATA U( 16),U( 18) /1.060460290274525D+01, 1.280182748008147D+01/
      DATA U( 20),U( 22) /1.510441257307552D+01, 1.750230784587389D+01/
      DATA U( 24),U( 26) /1.998721449566189D+01, 2.255216385312342D+01/
      DATA U( 28),U( 30) /2.519122118273868D+01, 2.789927138384089D+01/
      DATA U( 32),U( 34) /3.067186010608067D+01, 3.350507345013689D+01/
      DATA U( 36),U( 38) /3.639544520803305D+01, 3.933988418719949D+01/
      DATA U( 40),U( 42) /4.233561646075349D+01, 4.538013889847691D+01/
      DATA U( 44),U( 46) /4.847118135183522D+01, 5.160667556776437D+01/
      DATA U( 48),U( 50) /5.478472939811232D+01, 5.800360522298052D+01/
      DATA U( 52),U( 54) /6.126170176100200D+01, 6.455753862700633D+01/
      DATA U( 56),U( 58) /6.788974313718153D+01, 7.125703896716801D+01/
      DATA U( 60),U( 62) /7.465823634883016D+01, 7.809222355331531D+01/
      DATA U( 64),U( 66) /8.155795945611504D+01, 8.505446701758152D+01/
      DATA U( 68),U( 70) /8.858082754219768D+01, 9.213617560368709D+01/
      DATA U( 72),U( 74) /9.571969454214320D+01, 9.933061245478743D+01/
      DATA U( 76),U( 78) /1.029681986145138D+02, 1.066317602606435D+02/
      DATA U( 80),U( 82) /1.103206397147574D+02, 1.140342117814617D+02/
      DATA U( 84),U( 86) /1.177718813997451D+02, 1.215330815154386D+02/
      DATA U( 88),U( 90) /1.253172711493569D+02, 1.291239336391272D+02/
      DATA U( 92),U( 94) /1.329525750356163D+02, 1.368027226373264D+02/
      DATA U( 96),U( 98) /1.406739236482343D+02, 1.445657439463449D+02/
      DATA U(100),U(102) /1.484777669517730D+02, 1.524095925844974D+02/
      DATA U(104),U(106) /1.563608363030788D+02, 1.603311282166309D+02/
      DATA U(108),U(110) /1.643201122631952D+02, 1.683274454484277D+02/
      DATA U(112),U(114) /1.723527971391628D+02, 1.763958484069974D+02/
      DATA U(116),U(118) /1.804562914175438D+02, 1.845338288614495D+02/
      DATA U(120),U(122) /1.886281734236716D+02, 1.927390472878449D+02/
      DATA U(124),U(126) /1.968661816728900D+02, 2.010093163992815D+02/
      DATA U(128),U(130) /2.051681994826412D+02, 2.093425867525368D+02/
      DATA U(132),U(134) /2.135322414945633D+02, 2.177369341139542D+02/
      DATA U(136),U(138) /2.219564418191303D+02, 2.261905483237276D+02/
      DATA U(140),U(142) /2.304390435657770D+02, 2.347017234428183D+02/
      DATA U(144),U(146) /2.389783895618343D+02, 2.432688490029827D+02/
      DATA U(148),U(150) /2.475729140961869D+02, 2.518904022097232D+02/
      DATA U(152),U(154) /2.562211355500095D+02, 2.605649409718632D+02/
      DATA U(156),U(158) /2.649216497985528D+02, 2.692910976510198D+02/
      DATA U(160),U(162) /2.736731242856937D+02, 2.780675734403661D+02/
      DATA U(164),U(166) /2.824742926876304D+02, 2.868931332954270D+02/
      DATA U(168),U(170) /2.913239500942703D+02, 2.957666013507606D+02/
      DATA U(172),U(174) /3.002209486470141D+02, 3.046868567656687D+02/
      DATA U(176),U(178) /3.091641935801469D+02, 3.136528299498791D+02/
      DATA U(180),U(182) /3.181526396202093D+02, 3.226634991267262D+02/
      DATA U(184),U(186) /3.271852877037752D+02, 3.317178871969285D+02/
      DATA U(188),U(190) /3.362611819791985D+02, 3.408150588707990D+02/
      DATA U(192),U(194) /3.453794070622669D+02, 3.499541180407702D+02/
      DATA U(196),U(198) /3.545390855194408D+02, 3.591342053695754D+02/
      DATA U(200),U(202) /3.637393755555635D+02, 3.683544960724047D+02/
 
      LCG=.FALSE.
      GO TO 7
 
      ENTRY DCLEBG(A1,B1,C1,X1,Y1,Z1)
      LCG=.TRUE.
 
    7 H=0
      IA=NINT(2*A1)
      IB=NINT(2*B1)
      IC=NINT(2*C1)
      IX=NINT(2*X1)
      IY=NINT(2*Y1)
      IZ=NINT(2*Z1)
      IF(IA .LT. 0 .OR. IB .LT. 0 .OR. IC .LT. 0) GO TO 99
      IF(MOD(IA+IB+IC,2) .NE. 0) GO TO 99
      JX=ABS(IX)
      JY=ABS(IY)
      JZ=ABS(IZ)
      IF(IA .LT. JX .OR. IB .LT. JY .OR. IC .LT. JZ) GO TO 99
      IF(MOD(IA+JX,2) .NE. 0 .OR. MOD(IB+JY,2) .NE. 0) GOTO 99
      IF(MOD(IC+JZ,2) .NE. 0) GO TO 99
      IF(LCG) THEN
       IZ=-IZ
       J0=0
       F=SQRT((IC+1)*R1)
      ELSE
       J0=IA-IB-IZ
       F=1
      ENDIF
      IF(IX+IY+IZ .NE. 0 .OR. MOD(J0,2) .NE. 0) GO TO 99
      K0=IA+IB+IC+2
      K1=IA+IB-IC
      K2=IA-IB+IC
      K3=IB+IC-IA
      IF(K1 .LT. 0 .OR. K2 .LT. 0 .OR. K3 .LT. 0) GO TO 99
      K4=IA+IX
      K5=IB+IY
      K6=IC+IZ
      K7=IA-IX
      K8=IB-IY
      K9=IC-IZ
      K10=IB-IC-IX
      K11=IA-IC+IY
      KA=MAX(0,K10,K11)
      KZ=MIN(K1,K5,K7)
      W=HF*(U(K1)+U(K2)+U(K3)+U(K4)+U(K5)+U(K6)+U(K7)+U(K8)+U(K9)-U(K0))
      S=0
      Q=(-1)**((KA+J0)/2)
      DO 1 K = KA,KZ,2
      S=S+Q*EXP(W-(U(K)+U(K1-K)+U(K5-K)+U(K7-K)+U(K-K10)+U(K-K11)))
    1 Q=-Q
      H=F*S
      GO TO 99
 
      ENTRY DWIG6J(A1,B1,C1,X1,Y1,Z1)
 
      LJN=.FALSE.
      LRC=.FALSE.
      A=A1
      B=B1
      C=C1
      X=X1
      Y=Y1
      Z=Z1
      GO TO 9
 
      ENTRY DRACAW(A1,B1,C1,X1,Y1,Z1)
 
      LJN=.FALSE.
      LRC=.TRUE.
      GO TO 8
 
      ENTRY DJAHNU(A1,B1,C1,X1,Y1,Z1)
 
      LJN=.TRUE.
      LRC=.FALSE.
    8 A=A1
      B=B1
      C=Y1
      X=X1
      Y=C1
      Z=Z1
 
    9 H=0
      IA=NINT(2*A)
      IB=NINT(2*B)
      IC=NINT(2*C)
      IF(IA .LT. 0 .OR. IB .LT. 0 .OR. IC .LT. 0) GO TO 99
      IX=NINT(2*X)
      IY=NINT(2*Y)
      IZ=NINT(2*Z)
      IF(IX .LT. 0 .OR. IY .LT. 0 .OR. IZ .LT. 0) GO TO 99
      IABC=IA+IB+IC
      IAYZ=IA+IY+IZ
      IF(MOD(IABC,2) .NE. 0 .OR. MOD(IAYZ,2) .NE. 0) GOTO 99
      IXBZ=IX+IB+IZ
      IXYC=IX+IY+IC
      IF(MOD(IXBZ,2) .NE. 0 .OR. MOD(IXYC,2) .NE. 0) GOTO 99
      K1=IA+IB-IC
      K2=IA-IB+IC
      K3=IB+IC-IA
      IF(K1 .LT. 0 .OR. K2 .LT. 0 .OR. K3 .LT. 0) GO TO 99
      K4=IA+IY-IZ
      K5=IA-IY+IZ
      K6=IY+IZ-IA
      IF(K4 .LT. 0 .OR. K5 .LT. 0 .OR. K6 .LT. 0) GO TO 99
      K7=IX+IB-IZ
      K8=IX-IB+IZ
      K9=IB+IZ-IX
      IF(K7 .LT. 0 .OR. K8 .LT. 0 .OR. K9 .LT. 0) GO TO 99
      K10=IX+IY-IC
      K11=IX-IY+IC
      K12=IY+IC-IX
      IF(K10 .LT. 0 .OR. K11 .LT. 0 .OR. K12 .LT. 0) GO TO 99
      IABXY=IA+IB+IX+IY
      IBCYZ=IB+IC+IY+IZ
      ICAZX=IC+IA+IZ+IX
      KA=MAX(IABC,IAYZ,IXBZ,IXYC)
      KZ=MIN(IABXY,IBCYZ,ICAZX)
      J1=KA
      IF(LRC .OR. LJN) J1=KA+IABXY
      W=HF*(U(K1)+U(K2)+U(K3)-U(IABC+2)+U(K4)+U(K5)+U(K6)-U(IAYZ+2)+
     1      U(K7)+U(K8)+U(K9)-U(IXBZ+2)+U(K10)+U(K11)+U(K12)-U(IXYC+2))
      S=0
      Q=(-1)**(J1/2)
      DO 2 K = KA,KZ,2
      S=S+Q*EXP(W+U(K+2)-(U(K-IABC)+U(K-IAYZ)+U(K-IXBZ)+U(K-IXYC)+
     1            U(IABXY-K)+U(IBCYZ-K)+U(ICAZX-K)))
    2 Q=-Q
      H=S
      IF(LJN) H=SQRT(((IC+1)*(IZ+1))*R1)*H
 
   99 DWIG3J=H
      RETURN
      END