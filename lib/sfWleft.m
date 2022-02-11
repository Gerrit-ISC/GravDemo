function [sys,x0,str,ts]=sfWleft(t,x,u,flag,wn,wperc0,campl)

persistent W CA U

switch flag,
    case 0 %init
        W.n0=wn;
        W.perc0=wperc0;
        CA.mpl=campl;
        U=0;
        str = [];
        x0=[0];
        sys = [0 1 2 1 0 0 1];%xc xd y u 0 0 1
        ts = 0.001;
    
	case 2 %derivative
        U=u;
        sys=[0];
    
	case 3 %output
        
		n=0;
        mass=0;
            
        theta=U/360;
        
        n=W.n0;
        dist=theta/10.99206;
        I=round(dist);
        D=dist-I;
        n=n+I;
        DeltaD=W.perc0+D*100;
        
        if DeltaD>100
            n=n+1;
            DeltaD=DeltaD-100;
        elseif DeltaD<0
            n=n-1;
            DeltaD=100+DeltaD;
        end
        
        l = [20720;20983;21247;21511;21775]'*0.001;
        
        if n-1>0
            for i=1:1:n-1
                mass=mass+l(i)*CA.mpl;
            end
        end
        
        if n>0
            mass = mass + l(n)*DeltaD*CA.mpl/100;
        end
         if isempty(n)
            n=0;
        end
       
        sys=[mass;n];
    
	otherwise
        sys=[];
        
end