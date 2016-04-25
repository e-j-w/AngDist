	PROGRAM ang_dist

c	Program calculates directional distributions of gamma rays
c	according to eq. 12.197 of W.D Hamilton et al., 'The Electromagnetic
c	Interaction in Nuclear Spectroscopy'

c	Statistical tensor and F-coefficents are calculated using existing 
c	code by K. Starosta.
	
c	Declare functions
	REAL*8 B,F1,A,LPI
c	Declare variables
	REAL*8 lambda,Ival,sigmaj
	REAL*8 l,lprime,I0,I1,delta
	REAL*8 lval,hval,dist_val
	
	
	WRITE(*,1)''
	WRITE(*,1)'ANGULAR DISTRIBUTION CALCULATOR'
	WRITE(*,1)'-------------------------------'
	WRITE(*,1)''
1	FORMAT('',A31)

	WRITE(*,2)'Enter [lambda, I, sigma/I]:'
2	FORMAT('',A27)

	read(*,*)lambda,Ival,sigmaj

c	lambda=4.
c	Ival=4.
c	sigmaj=0.2

	
c	WRITE(*,100)B(lambda,Ival,sigmaj)
c100	FORMAT('B = ',F15.10)

	WRITE(*,3)'Enter [lambda, L, L`, I0, I1]:'
3	FORMAT('',A30)

	read(*,*)lambda,l,lprime,I0,I1
c	lambda=4.
c	l=2.
c	lprime=2.
c	I0=4.
c	I1=4.
c	For now, assume no mixing
	delta=0.0

c	WRITE(*,101)F1(lambda,l,lprime,I0,I1)
c101	FORMAT('F1 = ',F15.10)

	WRITE(*,*)"B = ",B(lambda,Ival,sigmaj)
	WRITE(*,*)"F = ",F1(lambda,l,lprime,I0,I1)
	WRITE(*,*)"A = ",A(lambda,l,lprime,I0,I1,delta)
	
	WRITE(*,4)'Enter angular range in deg [low,high]:'
4	FORMAT('',A38)

	read(*,*)lval,hval
	
	

	dist_val=0.
	do i=0,lambda,2
		WRITE(*,*)"i = ",i
c		Need to execute this on its own line rather than a write line 
c		to avoid recursive write statements (since there are write 
c		statements in this function)
		dist_val=dist_val+LPI(i,cos(lval*3.14159265359/180),cos(hval*3.14159265359/180))*B(lambda,Ival,sigmaj)*A(lambda,l,lprime,I0,I1,delta)
	end do

	WRITE(*,*)"Angular distribution coefficient = ",dist_val

	END
	
	
c	Function calculates angular distribution coefficients
c	according to eq. 12.185, pg. 542
	REAL*8 FUNCTION A(LB,L,L1,IF,II,delta)

	IMPLICIT NONE

C	Declare functions
	REAL*8 F1

C	Declare global variables
	REAL*8 LB,L,L1,IF,II,delta

	A=(F1(LB,L,L1,IF,II) + 2.*delta*F1(LB,L,L1+1,IF,II) + delta*delta*F1(LB,L+1,L1+1,IF,II))/(1 + delta*delta)
	
	RETURN
	END

	
c	Function calculates integral of the Legendre polynomial of the
c	specified order over the range low->high
	REAL*8 FUNCTION LPI(order,low,high)

	IMPLICIT NONE

C	Declare global variables
	REAL*8 low,high
	integer order
	REAL*8 tmp
	
	LPI=0.

c	Swap low and high if needed	
c	if(low.gt.high) then
c		tmp=low
c		low=high
c		high=tmp
c	endif
	
	if(order.eq.0) then
		LPI=high-low
		RETURN
	else if(order.eq.1) then
		LPI=((high**2)/2) - ((low**2)/2)
		RETURN
	else if(order.eq.2) then
		LPI=0.5*((high**3 - high) - (low**3 - low))
		RETURN
	else if(order.eq.3) then
		LPI=0.5*(((5/4)*high**4 - (3/2)*high**2) - ((5/4)*low**4 - (3/2)*low**2))
		RETURN
	else if(order.eq.4) then
		LPI=(0.125)*((7*high**5 - 10*high**3 + 3*high) - (7*low**5 - 10*low**3 + 3*low))
		RETURN
	else
		WRITE(*,*)"Only orders up to 4 are implemented!"
	endif
	
	
	RETURN
	END
