import math 

Rho=1
Nu=0.804*10e-6
form='normal_section'

def wetted_surface (L, T, B, Cm, Cb, Cwp):
    a=L*((2*T+B)*pow(Cm, 0.5))
    b=(0.453+(0.4425*Cb)-(0.2862*Cm)-(0.003467*(B/T))+(0.3696*Cwp))#+(2.38*Abt/Cb)
    return (a*b)

# print(wetted_surface(6.049, 0.38, 1.13, 0.98, 0.81, 0.89))

def form_factor(L, LCB, Cp, T, form, B ):
    Lr=(1-Cp+(0.06*Cp*LCB)/(4*Cp-1))*L;
    C12=0;
    k=T/L;
    if k>=0.05:
        C12=pow(k,0.2228446)
    elif k>=0.02 and k<0.05:
        C12=(48.20*(pow(k-0.02, 2.078))+0.479948)
    else:
        C12=0.479948
    dict1={'v_shaped':-10,'normal_section':0,'u_shaped':10}
    for form in dict1.keys():
        Cstern=dict1.get(form);
        break
    C13=1+(0.003*Cstern)
    return (C13*(0.93+(C12*(pow(B/Lr, 0.92497))*(pow(0.95-Cp, -0.521448))*(pow(1-Cp+(0.0225*LCB), 0.6903)))))
    
# print(form_factor(6.049, 0.0299, 0.83, 0.38, 'normal_section', 1.13))  

def fr_res_coeff(V, L, Nu):
    Re=reynolds_no(V, L, Nu)
    a=(math.log(Re,10))
    Cf=0.075/(pow((a-2),2))
    return Cf 

def reynolds_no(V, L, Nu):
    return V*L/Nu

def froude_no(V,L):
    return V/pow(L*9.81, 0.5)

# print ("froude number is" ,froude_no(2.57, 6.049))
    
    
def frictional_resistance( V, L, T, B, Cm, Cb, Cwp):
    S=wetted_surface(L, T, B, Cm, Cb, Cwp)
    Cf=fr_res_coeff(V, L, Nu)
    return 0.5*Rho*V*V*S*Cf;
# print(frictional_resistance( 2.57, 6.049, 0.38, 1.13, 0.98, 0.81, 0.89))

# def appendage_resistance(Rho,K,V,L,Nu,T, B, Cm, Cb, Cwp, Abt):
#
#     Cf=fr_res_coeff(V, L, Nu)
#     S=wetted_surface(L, T, B, Cm, Cb, Cwp, Abt)
#     return 0.5*Rho*Cf*S*(1+K);

def wave_resistance(V, B, L, T, Cwp, Cp, lcb, Lr, disp):
    a=B/L
    if(a<0.11):
        C7=0.229577*pow(a, 0.33333)
    elif(a>=0.11 and a<=0.25):
        C7=a
    else:
        C7=0.5-(0.0625*a)
    e=2.71828
    b=pow(a, -0.80856)*pow(1-Cwp, 0.30484)*pow(1-Cp-(0.0225*lcb), 0.6367)*pow(Lr/B, 0.34574)*pow(100*disp/(L*L*L), 0.16302)
    iE=1+pow(89, -b)
    C1=2223105*pow(C7, 3.78613)*pow(T/B, 1.07961)*pow((90-iE), -1.37565)
    
    if(L/B<12):
        Lam=1.446*Cp-(0.03*L/B)
    else:
        Lam=1.446*Cp-0.36
    if(Cp>=0.80):
        C16=1.73014-(0.7067*Cp)
    else:
        C16=8.07981*Cp-(13.8673*Cp*Cp)+(6.984388*Cp*Cp*Cp)
    
    
    m1=(0.0140407*L/T)-(1.75254*pow(disp, 1/3)/L)-(4.79323*B/L)-C16
    c=L*L*L/disp
    if(c<512):
        C15=-1.69385
    elif(c>=512 and c<1727):
        C15=-1.69385+(((L/pow(disp, 1/3))-8.0)/2.36)
    
    Fn=froude_no(V, L)
    m2=C15*Cp*Cp*pow(e, -0.1/(Fn*Fn))
    d=-0.9
    return C1*disp*Rho*9.81*pow(e,  (m1*pow(Fn, d))+(m2*math.cos(Lam/(Fn*Fn))));

# print(wave_resistance(2.57,1.13, 6.049, 0.38, 0.89, 0.83, 2.99, 1.032, 2.090))


# def wind_resistance():
#     return 11.36;

def calc_Cb(disp,length,breadth,draft):
        return disp/(length*breadth*draft)

class resistance:
    def __init__(self,length,breadth,draft,V,Cp,Cm,Disp,Cwp,LCB,form):
        self.length=length;
        self.breadth=breadth;
        self.draft=draft
        self.V=V                        # speed in m/s
        self.Cp=Cp                      # prismatic coefficient
        self.Cm=Cm                      # Mid-ship Coefficient                      
        self.disp=Disp                  # Displacement in ton
        self.Cwp=Cwp                    # Waterplane Coefficient
        self.LCB=LCB                    # Longitudinal Centre of Buoyancy
        self.form=form
        self.Cb=calc_Cb(self.disp,self.length,self.breadth,self.draft)
                                        # Block Coefficient
                                        
    def total_resistance(self):
        Rf=frictional_resistance( self.V, self.length, self.draft, self.breadth, self.Cm, self.Cb, self.Cwp)
        Lr=(1-self.Cp+(0.06*self.Cp*self.LCB)/(4*self.Cp-1))*self.length;
        Rw=wave_resistance(self.V, self.breadth, self.length, self.draft, self.Cwp, self.Cp, self.LCB, Lr, self.disp)
        
        Rt=form_factor(self.length, self.LCB, self.Cp, self.draft, self.form, self.breadth)*Rf+Rw
        
        print("Total Resistance = ",Rt,"KN")
        print("Total Resistance (kgf) = ",Rt*102,"kgf")
        print("Effective Power = ",Rt*self.V,"KW")
        
        
ship1=resistance(6.049, 1.13, 0.38, 2.57, 0.83, 0.98, 2.09, 0.89, 0.004,'normal_section')
ship1.total_resistance()
ship2=resistance(4.20, 1.65, 0.97,  1.02, 0.84, 0.98, .1, 309.16/386.97, 2.99,'normal_section')
ship2.total_resistance()
