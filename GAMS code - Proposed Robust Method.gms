sets
a        /1*2/
i       /1*8/
j       /1*2/
k       /1*2/
l       /1*2/
r       /1*2/
w       /1*2/
;

scalars
b        /8/
alpha    /0.18/
beta     /0.69/
M        /10000000000/
teta     /0.5/
mu       /0.3/
;


parameters
gamma(r) /1 0.04,2 0.04/

wh(w)    /1 2686152,2 1224052/

cw(w)    /1 667585,2 984853/

qrij(r)  /1 10.35,2 13.85/

qrik(r)  /1 10.80,2 12.30 /

qrkl(r)  /1 25.48,2 25.53/

qrjw(r)  /1 18.30,2 19.25/

qrlw(r)  /1 38.79,2 36.31/

Crik(r)  /1 87,2 65 /

Crij(r)  /1 51,2 73 /

Crkl(r)  /1 139,2 144/

Crlw(r)  /1 81,2 72 /

Crjw(r)  /1 196,2 177/

nn(a)    /1 3743,2   10708/
;

parameter dij(i,j);

$ call gdxxrw Samples.xlsx par dij rng=Dij!B2:C9 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load dij
$gdxin

parameter dik(i,k);

$ call gdxxrw Samples.xlsx par dik rng=Dik!B2:C9 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load dik
$gdxin

parameter dkl(k,l);

$ call gdxxrw Samples.xlsx par dkl rng=Dkl!B2:C3 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load dkl
$gdxin

parameter dlw(l,w);

$ call gdxxrw Samples.xlsx par dlw rng=Dlw!B2:C3 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load dlw
$gdxin

parameter djw(j,w);

$ call gdxxrw Samples.xlsx par djw rng=Djw!B2:C3 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load djw
$gdxin

parameter daj(a,j);

$ call gdxxrw Samples.xlsx par daj rng=Daj!B2:C3 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load daj
$gdxin

parameter dal(a,l);

$ call gdxxrw Samples.xlsx par dal rng=Dij!B2:C3 rdim=1 cdim=1
$ gdxin Samples.gdx
$ load dal
$gdxin

parameters
po
pm
pp
fo(i)
fm(i)
fp(i)
so
sm
sp
nao(a)
nam(a)
nap(a)

;

po=normal(17,1);
pm=normal(23,1);
pp=normal(29,1);

fp(i)= normal(197,1);
fm(i)= normal(202,1);
fo(i)= normal(207,1);

sp= normal(0.2,1);
sm= normal(0.4,1);
so= normal(0.6,1);

nap(a)= normal(300,1);
nam(a)= normal(800,1);
nao(a)= normal(1300,1);



parameters
Bk
Bj
Bl
Bw
ej(j)
el(l)
gj(j)
gl(l)
Hj
Hl
ok
oj
ol
ow
ck
cj
cl
;

Bk=1025;
Bj=1650 ;
Bl=16500;
Bw=44000 ;
ej(j)=9 ;
el(l)=9;
gj(j)=1650;
gl(l)=44000;
Hj=0.22;
Hl=0.22;
ok=2;
oj=45.5;
ol=45.3;
ow=2;
ck=1;
cj=165000;
cl=802000;

binary variables
zk(k)
zj(j)
zl(l)
zw(w)
;

free variables
z1
z2
z3
zt
;
integer variables
nrij(r,i,j)
nrik(r,i,k)
nrkl(r,k,l)
nrjw(r,j,w)
nrlw(r,l,w)
paj(a,j)
pal(a,l)
pj(j)
pl(l)
;

positive variables
xij(i,j)
xik(i,k)
xkl(k,l)
yjw(j,w)
ylw(l,w)
yw(w);

equations
obj1
obj2
obj3
TotalObj
co1
co2
co3
co4
co5
co6
co7
co8
co9
co10
co11
co12
co13
co14
*co15
co16
co17
co18
co19
co20
co21
co22
co23
co24
co25
co26
co27
co28
co29
co30
co31
co32
co33
co34
co35

;


obj1     ..        z1 =e= [(po+4*pm+pp)/6]*(sum((i,j),xij(i,j))+sum((i,k),xik(i,k)))+ mu*[(pp-po)/4.9]*(sum((i,j),xij(i,j))+sum((i,k),xik(i,k)))
                         +sum((r,i,j),Crij(r)*dij(i,j)*nrij(r,i,j))+sum((r,i,k),Crik(r)*dik(i,k)*nrik(r,i,k))
                         +sum((r,k,l),Crkl(r)*dkl(k,l)*nrkl(r,k,l))+sum((r,j,w),Crjw(r)*djw(j,w)*nrjw(r,j,w))+sum((r,l,w),Crlw(r)*dlw(l,w)* nrlw(r,l,w))
                         +sum(k,ck*zk(k))+sum(l,cl*zl(l))+ sum(j,cj*zj(j))+ sum(w,cw(w)*zw(w))
                         +[ok*sum((i,k),xik(i,k))]+[oj*sum((j,w),yjw(j,w))]
                         +[ol*sum((l,w),ylw(l,w))]+ [ow*sum(w,yw(w))];

obj2     ..      z2 =e= alpha*[(sum((a,j),daj(a,j)*paj(a,j))+sum((a,l),dal(a,l)*pal(a,l)))]
                        +beta*[sum(w,yw(w))]
                               +sum((r,i,j),gamma(r)*dij(i,j)*nrij(r,i,j))+sum((r,i,k),gamma(r)*dik(i,k)*nrik(r,i,k))
                               +sum((r,k,l),gamma(r)*dkl(k,l)*nrkl(r,k,l))+sum((r,j,w),gamma(r)*djw(j,w)*nrjw(r,j,w))+sum((r,l,w),gamma(r)*dlw(l,w)*nrlw(r,l,w)) ;

obj3     ..      z3 =e=  sum((j,w),(hj*yjw(j,w))/b)+sum((l,w),(hl*ylw(l,w))/b)+[9*(sum(j,zj(j))+sum(l,zl(l)))];



TotalObj        ..      ZT =e= 1/3*[z1/927428]+1/3*[z2/802]-1/3*[z3/193];



*co1  ..      sum(j,xij(i,j))+  sum(k,xik(i,k)) =g= f;
co1(i)     ..      sum(j,xij(i,j))+  sum(k,xik(i,k)) =g= teta*(2/3*fm(i)+1/3*fo(i))+(1-teta)*(2/3*fm(i)+1/3*fp(i));

co2(k)  ..      sum(i,xik(i,k)) =l= Bk;
co3(j)  ..      sum(i,xij(i,j)) =l= Bj;
co4(l)  ..      sum(k,xkl(k,l)) =l= Bl;
co5(w) ..      sum(l,ylw(l,w))+sum(j,yjw(j,w)) =l= Bw;

co6(k)  ..      sum(i,xik(i,k)) =E= sum(l,xkl(k,l));
co7(w)  ..      sum(j,yjw(j,w))+ sum(l,ylw(l,w)) =e= yw(w);

*co8(j)  ..    s*sum(i,xij(i,j)) =e=  sum(w,yjw(j,w));
co8(j)   ..      [(1-teta/2)*(2/3*sm+1/3*so)+(teta/2)*(2/3*sm+1/3*sp)]*sum(i,xij(i,j)) =g=  sum(w,yjw(j,w));
co9(j)   ..      [(teta/2)*(2/3*sm+1/3*so)+(1-teta/2)*(2/3*sm+1/3*sp)]*sum(i,xij(i,j)) =l=  sum(w,yjw(j,w));


*co10(l)  ..    s*sum(k,xkl(k,l)) =e=  sum(w,ylw(l,w));
co10(l)   ..    [(1-teta/2)*(2/3*sm+1/3*so)+(teta/2)*(2/3*sm+1/3*sp)]* sum(k,xkl(k,l)) =g=  sum(w,ylw(l,w));
co11(l)   ..    [(teta/2)*(2/3*sm+1/3*so)+(1-teta/2)*(2/3*sm+1/3*sp)]* sum(k,xkl(k,l)) =l=  sum(w,ylw(l,w));



co14(a) ..   -1*sum(j,paj(a,j))+sum(l,pal(a,l)) =g= -1*[teta*(2/3*nam(a)+1/3*nao(a))+(1-teta)*(2/3*nam(a)+1/3*nap(a))];
*co14(a) ..      sum(j,paj(a,j)) =L= nn(a);
*co15(a) ..      sum(l,pal(a,l)) =L= nn(a);

co12(j) ..      sum(w,yjw(j,w)) =l= gj(j);
co13(l) ..      sum(w,ylw(l,w)) =l= gl(l);
co16(r,i,j)      ..      nrij(r,i,j) =g= xij(i,j)/qrij(r);
co17(r,i,k)      ..      nrik(r,i,k) =g= xik(i,k)/qrik(r) ;
co18(r,k,l)      ..      nrkl(r,k,l) =g= xkl(k,l)/qrkl(r) ;
co19(r,j,w)      ..      nrjw(r,j,w) =g= yjw(j,w)/qrjw(r) ;
co20(r,l,w)      ..      nrlw(r,l,w) =g= ylw(l,w)/qrlw(r) ;
co21(j)  ..    sum(a,paj(a,j))=e=ej(j)*zj(j)+[sum(w,hj*yjw(j,w))/b];
co22(l)  ..     sum(a,pal(a,l))=e=el(l)*zl(l)+[sum(w,hl*ylw(l,w))/b];
co23(j)  ..      sum(a,paj(a,j)) =l= M*zj(j);
co24(l)  ..      sum(a,pal(a,l)) =l= M*zl(l);
co25(k) ..      sum(i,xik(i,k)) =l= M*zk(k);
co26(j) ..      sum(i,xij(i,j)) =l= M*zj(j);
co27(l) ..      sum(k,xkl(k,l)) =l= M*zl(l);
co28(w) ..      sum(j,yjw(j,w)) =l= M*zw(w);
co29(w) ..      sum(l,ylw(l,w)) =l= M*zw(w);
co30(w) ..      sum(j,yjw(j,w)) =l= wh(w);
co31(w) ..      sum(l,ylw(l,w)) =l= wh(w);
co32(k)  ..     zk(k) =l= 1;
co33(j)  ..     zj(j) =l= 1;
co34(l)  ..     zl(l) =l= 1;
co35(w) ..      zw(w) =l= 1;




model   PessimisticFuzzyRobust    /all/

options
optca=0
optcr=0
reslim=500000
limrow=10000
limcol=10000
MIP=CPLEX
decimals=5
;


solve PessimisticFuzzyRobust using MIP minimizing zt ;

display ZT.l,z1.l,z2.l,z3.l,zk.l,zj.l,zl.l,zw.l,paj.l,pal.l,xij.l,xik.l,xkl.l,yjw.l,ylw.l,yw.l,nrik.l,nrij.l,nrkl.l,nrlw.l,nrjw.l


