import sys
sys.path.append('./toolbox/')

import Deg2Solver as D2S

def LinearShooting(p,q,r,alpha,beta,a,b,h,n):
    """ 
    This functions solves BVP, and returns the function value and its derivative at all the mesh
    points , using the linear shooting and RK4 method.
    PROBLEM:
    Given: ypp(t,y,yp)=p(t)*yp(t)+q(t)*y(t)+r(t) ;
           y(a)=alpha ;
           y(b)=beta ;
           a+h*n=b
    """
    def f1(X,Y,YP):
        return (p(X)*YP)+(q(X)*Y)+r(X)
    
    def f2(X,Y,YP):
        return (p(X)*YP)+(q(X)*Y)
    
    (y1is,y1pis)=D2S.Deg2RK4(f1,a,alpha,0,h,n)
    (y2is,y2pis)=D2S.Deg2RK4(f2,a,0,1,h,n)
    yps=[]
    ys=[]
    for i in range(n+1):
        yps.append(y1pis[i] + (((beta - y1is[n]) / y2is[n]) * y2pis[i]))
    for i in range(n+1):
        ys.append(y1is[i] + (((beta - y1is[n]) / y2is[n]) * y2is[i]))
    return (ys,yps)

def shoot_nonlinear(a,b,alpha, beta, n, tol, M,f,fy,fyp):

    w1 = [0 for i in range(n+1)]
    w2 = [0 for i in range(n+1)]
    h = (b-a)/n
    k = 1
    TK = (beta - alpha)/(b - a)

    while k <= M:

        w1[0] = alpha
        w2[0] = TK
        u1    = 0
        u2    = 1

        for i in range(1,n+1):
            x = a + (i-1)*h     #step 5

            t = x + 0.5*(h)

            k11 = h*w2[i-1]     #step 6

            k12 = h*f(x,w1[i-1],w2[i-1])
            k21 = h*(w2[i-1] + (1/2)*k12)
            k22 = h*f(t, w1[i-1] + (1/2)*k11, w2[i-1] + (1/2)*k12)
            k31 = h*(w2[i-1] + (1/2)*k22)
            k32 = h*f(t, w1[i-1] + (1/2)*k21, w2[i-1] + (1/2)*k22)
            t   = x + h
            k41 = h*(w2[i-1]+k32)
            k42 = h*f(t, w1[i-1] + k31, w2[i-1] + k32)
            w1[i] = w1[i-1] + (k11 + 2*k21 + 2*k31 + k41)/6
            w2[i] = w2[i-1] + (k12 + 2*k22 + 2*k32 + k42)/6   
            kp11 = h*u2
            kp12 = h*(fy(x,w1[i-1],w2[i-1])*u1 + fyp(x,w1[i-1], w2[i-1])*u2)
            t    = x + 0.5*(h)
            kp21 = h*(u2 + (1/2)*kp12)
            kp22 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp11)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp12))
            kp31 = h*(u2 + (1/2)*kp22)
            kp32 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp21)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp22))
            t    = x + h
            kp41 = h*(u2 + kp32)
            kp42 = h*(fy(t, w1[i-1], w2[i-1])*(u1+kp31) + fyp(x + h, w1[i-1], w2[i-1])*(u2 + kp32))
            u1 = u1 + (1/6)*(kp11 + 2*kp21 + 2*kp31 + kp41)
            u2 = u2 + (1/6)*(kp12 + 2*kp22 + 2*kp32 + kp42)


        r = abs(w1[n] - beta)
        #print(r)
        if r < tol:
                xs = [a + i*h for i in range(n+1)]
                return (xs,w1,w2)            

        TK = TK -(w1[n]-beta)/u1

        k = k+1


    print("Maximum number of iterations exceeded")   
    return
    
def ypp(X,Y,YP):
    return (1/8)(32+(2(X**3))-Y*YP)

def fy(X,Y,YP):
    return (-1/8)*YP

def fyp(X,Y,YP):
    return (-1/8)*Y

(xs,ys,ypps)=BVPS.shoot_nonlinear(1,3,17,(43/3),20,0.00001,100,ypp,fy,fyp)
print("Soln0=",ys)