import numpy as np

class Constant_Volume():
    def __init__(self,T_1,T_2,P_1,P_2,q_12,c_v,s_1,s_2):
        self.T_1 = T_1
        self.T_2 = T_2
        self.P_1 = P_1
        self.P_2 = P_2
        self.q_12 = q_12
        self.c_v = c_v
        self.s_1 = s_1
        self.s_2 = s_2

    def temp_and_pressure(self):
        # Checks if there are more than 1 unknown for the equation
        #check = [self.T_1, self.T_2, self.P_1, self.P_2]
        #null = 0
        #for a in check:
        #    if a == None:
        #        null += 1

        #if null < 1:
        #    #return [[self.T_1, self.T_2, self.P_1, self.P_2], 'Too many unknown values']


        #Checks which value is the unknown and calculates
        if self.T_1 == None:
            self.T_1 = self.T_2*self.P_1/self.P_2
            #return [[self.T_1, self.T_2, self.P_1, self.P_2], None]

        elif self.T_2 == None:
            self.T_2 = self.T_1*self.P_2/self.P_1
            #return [[self.T_1, self.T_2, self.P_1, self.P_2], None]

        elif self.P_1 == None:
            self.P_1 = self.P_2*self.T_1/self.T_2
            #return [[self.T_1, self.T_2, self.P_1, self.P_2], None]

        elif self.P_2 == None:
            self.P_2 = self.P_1*self.T_2/self.T_1
            #return [[self.T_1, self.T_2, self.P_1, self.P_2], None]



    def heat_transferred(self):
        # Checks if there are more than 1 unknown for the equation
        #check = [self.q_12, self.c_v, self.T_1, self.T_2]
        #null = 0
        #for a in check:
        #    if a == None:
        #        null += 1

        #if null < 1:
        #    return [[self.q_12, self.c_v, self.T_1, self.T_2], 'Too many unknown values']

        #Checks which value is the unknown and calculates
        if self.q_12 == None:
            self.q_12 = self.c_v*(self.T_2-self.T_1)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]

        elif self.c_v == None:
            self.c_v = self.q_12/(self.T_2-self.T_1)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]

        elif self.T_1 == None:
            self.T_1 = self.T_2 - self.q_12/self.c_v
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]

        elif self.T_2 == None:
            self.T_2 = self.T_1 + self.q_12/self.c_v
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]



       
    def entropy_change(self):
        # Checks if there are more than 1 unknown for the equation
        #check = [self.q_12, self.c_v, self.T_1, self.T_2, self.s_1, self.s_2]
        #null = 0
        #for a in check:
        #    if a == None:
        #        null += 1

        #if null < 1:
        #    return [[self.q_12, self.c_v, self.T_1, self.T_2, self.s_1, self.s_2], 'Too many unknown values']


        if self.s_1 == None:
            self.s_1 = self.s_2 - self.c_v*np.log(self.T_2/self.T_1)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]
        
        elif self.s_2 == None:
            self.s_2 = self.s_1 + self.c_v*np.log(self.T_2/self.T_1)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]

        elif self.c_v == None:
            self.c_v = (self.s_2 - self.s_1)/(np.log(self.T_2/self.T_1))
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]
        
        elif self.T_2 == None:
            self.T_2 = self.T_1*np.exp((self.s_2-self.s_1)/self.c_v)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]

        elif self.T_1 == None:
            self.T_1 = self.T_2/np.exp((self.s_2-self.s_1)/self.c_v)
            #return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2], None]


    def equation_finder(self):
        called = None
        for i in range(0,3):
            checka = [self.q_12, self.c_v, self.T_1, self.T_2, self.s_1, self.s_2]
            null_entropy = 0
            for a in checka:
                if a == None:
                    null_entropy += 1
            if null_entropy == 1:
                Constant_Volume.entropy_change(self)
                called = 1

            checkb = [self.T_1, self.T_2, self.P_1, self.P_2]
            null_temp_press = 0
            for a in checkb:
                if a == None:
                    null_temp_press += 1
                    

            if null_temp_press == 1:
                Constant_Volume.temp_and_pressure(self)
                called =1

            checkc = [self.q_12, self.c_v, self.T_1, self.T_2]
            null_heat = 0
            for a in checkc:
                if a == None:
                    null_heat += 1

            if null_heat == 1:
                Constant_Volume.heat_transferred(self)
                called =1


        if called == None:
            return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2],['Unable to compute values with variables']]
        else:
            return [[self.T_1,self.T_2,self.P_1,self.P_2,self.q_12,self.c_v,self.s_1,self.s_2],'N/A']
          


class Constant_Pressure():
    def __init__(self,T_1,T_2,v_1,v_2,q_12,w_12,c_v,c_p,R,s_1,s_2,p):
        self.T_1 = T_1
        self.T_2 = T_2
        self.v_1 = v_1
        self.v_2 = v_2
        self.q_12 = q_12
        self.w_12 = w_12
        self.c_v = c_v
        self.c_p = c_p
        self.R = R
        self.s_1 = s_1
        self.s_2 = s_2
        self.p = p

    def temp_and_volume(self):
        if self.T_1 == None:
            self.T_1 = self.T_2*self.v_1/self.v_2
            return self.T_1

        elif self.T_2 == None:
            self.T_2 = self.T_1*self.v_2/self.v_1
            return self.T_2

        elif self.v_1 == None:
            self.v_1 = self.v_2*self.T_1/self.T_2
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1*self.T_2/self.T_1
            return self.v_2


    def work_done_vol(self):
        if self.w_12 == None:
            self.w_12 = self.p*(self.v_2 - self.v_1)
            return self.w_12
        
        elif self.p == None:
            self.p = self.w_12/(self.v_2 - self.v_1)
            return self.p

        elif self.v_2 == None:
            self.v_2 = self.w_12/self.p + self.v_1
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = self.w_12/self.p - self.v_2
            return self.v_1

    def work_done_temp(self):
        if self.w_12 == None:
            self.w_12 = self.R*(self.T_2 - self.T_1)
            return self.w_12
        
        elif self.R == None:
            self.R = self.w_12/(self.T_2 - self.T_1)
            return self.R

        elif self.T_2 == None:
            self.T_2 = self.w_12/self.R + self.T_1
            return self.v_2

        elif self.T_1 == None:
            self.T_1 = self.w_12/self.R - self.T_2
            return self.T_1


    def heat_transfer_CvR(self):
        if self.q_12 == None:
            self.q_12 = self.c_v*(self.T_2 - self.T_1) + self.R*(self.T_2 - self.T_1)
            return self.q_12
        
        elif self.c_v == None:
            self.c_v = (self.q_12 - self.R*(self.T_2 - self.T_1))/(self.T_2 - self.T_1)
            return self.c_v
        
        elif self.R == None:
            self.R = (self.q_12 - self.c_v*(self.T_2 - self.T_1))/(self.T_2 - self.T_1)
            return self.R

        elif self.T_1 == None:
            self.T_1 = self.T_2 - self.q_12/(self.c_v+self.R)
            return self.T_1

        elif self.T_2 == None:
            self.T_2 = self.T_1 + self.q_12/(self.c_v+self.R)
            return self.T_2
        
    def heat_transfer_Cp(self):
        if self.q_12 == None:
            self.q_12 = self.c_p*(self.T_2 - self.T_1)
            return self.q_12
        
        elif self.c_p == None:
            self.c_v = self.q_12/(self.T_2 - self.T_1)
            return self.c_v
        
        elif self.T_1 == None:
            self.T_1 = self.T_2 - self.q_12/(self.c_p)
            return self.T_1

        elif self.T_2 == None:
            self.T_2 = self.T_1 + self.q_12/(self.c_p)
            return self.T_2


    def entropy_change(self):
        if self.s_1 == None:
            self.s_1 = self.s_2 - self.c_p*np.log(self.T_2/self.T_1)
            return self.s_1
    
        elif self.s_2 == None:
            self.s_2 = self.s_1 + self.c_p*np.log(self.T_2/self.T_1)
            return self.s_2

        elif self.c_p == None:
            self.c_p = (self.s_2 - self.s_1)/(np.log(self.T_2/self.T_1))
            return self.c_p
    
        elif self.T_2 == None:
            self.T_2 = self.T_1*np.exp((self.s_2-self.s_1)/self.c_p)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2/np.exp((self.s_2-self.s_1)/self.c_p)
            return self.T_1


    def equation_finder(self):
        called = None
        for i in range(0,6):
            checka = [self.T_1, self.T_2, self.v_1, self.v_2]
            null_temp_vol = 0
            for a in checka:
                if a == None:
                    null_temp_vol += 1
            if null_temp_vol == 1:
                Constant_Pressure.temp_and_volume(self)
                called = 1

            checkb = [self.w_12, self.p, self.v_1, self.v_2]
            null_work_done_vol = 0
            for a in checkb:
                if a == None:
                    null_work_done_vol += 1
                    

            if null_work_done_vol == 1:
                Constant_Pressure.work_done_vol(self)
                called =1

            checkc = [self.T_1, self.T_2, self.w_12, self.R]
            null_work_done_temp = 0
            for a in checkc:
                if a == None:
                    null_work_done_temp += 1

            if null_work_done_temp == 1:
                Constant_Pressure.work_done_temp(self)
                called =1

    
            checkd = [self.q_12, self.c_v, self.R, self.T_1, self.T_2]
            null_heat_CvR = 0
            for a in checkd:
                if a == None:
                    null_heat_CvR += 1

            if null_heat_CvR == 1:
                Constant_Pressure.heat_transfer_CvR(self)
                called =1

            checke = [self.q_12, self.c_p, self.T_1, self.T_2]
            null_heat_Cp = 0
            for a in checkd:
                if a == None:
                    null_heat_Cp += 1

            if null_heat_Cp == 1:
                Constant_Pressure.heat_transfer_Cp(self)
                called =1

            checkf = [self.s_1, self.s_2, self.c_p, self.T_2, self.T_1]
            null_entropy = 0
            for a in checkd:
                if a == None:
                    null_entropy += 1

            if null_entropy == 1:
                Constant_Pressure.entropy_change(self)
                called =1

                
        if called == None:
            return [[self.T_1,self.T_2,self.v_1,self.v_2,self.q_12,self.w_12,self.c_v,self.c_p,self.R,self.s_1,self.s_2,self.p],['Unable to compute values with variables']]
        else:
            return [[self.T_1,self.T_2,self.v_1,self.v_2,self.q_12,self.w_12,self.c_v,self.c_p,self.R,self.s_1,self.s_2,self.p],'N/A']



class Constant_Temperature():
    def __init__(self,T,v_1,v_2,P_1,P_2,q_12,w_12,c_v,c_p,R,s_1,s_2):
        self.T = T
        self.P_1 = P_1
        self.P_2 = P_2
        self.v_1 = v_1
        self.v_2 = v_2
        self.q_12 = q_12
        self.w_12 = w_12
        self.c_v = c_v
        self.c_p = c_p
        self.R = R
        self.s_1 = s_1
        self.s_2 = s_2

    def pressure_and_volume(self):
        if self.P_1 == None:
            self.P_1 = self.P_2*self.v_1/self.v_2
            return self.P_1

        elif self.P_2 == None:
            self.P_2 = self.P_1*self.v_2/self.v_1
            return self.P_2

        elif self.v_1 == None:
            self.v_1 = self.v_2*self.P_1/self.P_2
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1*self.P_2/self.P_1
            return self.v_2

    def entropy_change(self):
        if self.s_1 == None:
            self.s_1 = self.s_2 - self.R*np.log(self.v_2/self.v_1)
            return self.s_1

        elif self.s_2 == None:
            self.s_2 = self.s_1 + self.R*np.log(self.v_2/self.v_1)
            return self.s_2
       
        elif self.v_2 == None:
            self.v_2 = self.v_1*np.exp((self.s_2-self.s_1)/self.R)
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = self.v_2/np.exp((self.s_2-self.s_1)/self.R)
            return self.v_1


    def heat_transferred_entropy(self):
        if self.q_12 == None:
            self.q_12 = self.T*(self.s_2-self.s_1)
            return self.q_12

        elif self.T == None:
            self.T = self.q_12/(self.s_2-self.s_1)
            return self.T

        elif self.s_1 == None:
            self.s_1 = self.s_2 - self.q_12/self.T
            return self.s_1

        elif self.s_2 == None:
            self.s_2 = self.s_1 + self.q_12/self.T
            return self.s_2


    def heat_transferred_R(self):
        if self.q_12 == None:
            self.q_12 = self.R*self.T*np.log(self.v_2/self.v_1)
            return self.q_12

        elif self.T == None:
            self.T = self.q_12/(self.R*np.log(self.v_2/self.v_1))
            return self.T

        elif self.R == None:
            self.R = self.q_12/(self.T*np.log(self.v_2/self.v_1))

        elif self.v_2 == None:
            self.v_2 = self.v_1*np.exp((self.q_12)/self.R*self.T)
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = self.v_2/np.exp((self.q_12)/self.R*self.T)
            return self.v_1

    def work_done(self):
        if self.q_12 == None:
            self.q_12 = self.w_12
            return self.q_12
        
        elif self.w_12 == None:
            self.w_12 = self.q_12
            return self.w_12


    def equation_finder(self):
        called = None
        for i in range(0,5):
            checka = [self.P_1, self.P_2, self.v_1, self.v_2]
            null_pres_vol = 0
            for a in checka:
                if a == None:
                    null_pres_vol += 1
            if null_pres_vol == 1:
                Constant_Pressure.pressure_and_volume(self)
                called = 1

            checkb = [self.s_1, self.s_2, self.v_1, self.v_2]
            null_entropy = 0
            for a in checkb:
                if a == None:
                    null_entropy += 1
                    

            if null_entropy == 1:
                Constant_Pressure.entropy_change(self)
                called =1

            checkc = [self.q_12, self.T, self.s_1, self.s_2]
            null_heat_entropy = 0
            for a in checkc:
                if a == None:
                    null_heat_entropy += 1

            if null_heat_entropy == 1:
                Constant_Pressure.heat_transferred_entropy(self)
                called =1


            checkd = [self.q_12, self.T,self.R, self.v_1, self.v_2]
            null_heat_R = 0
            for a in checkd:
                if a == None:
                    null_heat_R += 1

            if null_heat_R == 1:
                Constant_Pressure.heat_transferred_R(self)
                called =1

            checke = [self.w_12, self.q_12]
            null_work = 0
            for a in checke:
                if a == None:
                    null_work += 1

            if null_work == 1:
                Constant_Pressure.work_done(self)
                called =1

                
        if called == None:
            return [[self.T,self.v_1,self.v_2,self.P_1,self.P_2,self.q_12,self.w_12,self.c_v,self.c_p,self.R,self.s_1,self.s_2],['Unable to compute values with variables']]
        else:
            return [[self.T,self.v_1,self.v_2,self.P_1,self.P_2,self.q_12,self.w_12,self.c_v,self.c_p,self.R,self.s_1,self.s_2],'N/A']




class Adiabatic():
    def __init__(self,T_1,T_2,v_1,v_2,P_1,P_2,w_12,c_v,c_p,R,s_1,s_2,gamma):
        self.T = T
        self.P_1 = P_1
        self.P_2 = P_2
        self.v_1 = v_1
        self.v_2 = v_2
        self.w_12 = w_12
        self.c_v = c_v
        self.c_p = c_p
        self.R = R
        self.s_1 = s_1
        self.s_2 = s_2
        self.gamma = gamma
        self.u_1 = u_1
        self.u_2 = u_2

    def reversible(self):
        if self.s_1 == None:
            self.s_1 = self.s_2
            return self.s_1

        elif self.s_2 == None:
            self.s_2 = self.s_1
            return self.s_2


    def gas_law_pres_vol(self):
        if self.P_2 == None:
            self.P_2 = self.P_1 * (self.v_1/self.v_2) ** self.gamma
            return self.P_2

        elif self.P_1 == None:
            self.P_1 = self.P_2 / ((self.v_1/self.v_2) ** self.gamma)
            return self.P_1

        elif self.v_1 == None:
            self.v_1 = self.v_2 * (self.P_1/self.P_2) ** self.gamma
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1 / ((self.P_1/self.P_2) ** ((self.gamma-1)/self.gamma))
            self.v_2

        elif self.gamma == None:
            self.gamma = np.log(self.P_2/self.P_1)/(np.log(self.v_1/self.v_2))
            return self.gamma

    def gas_law_temp_pres(self):
        if self.T_2 == None:
            self.T_2 = self.T_1 * (self.P_1/self.P_2) ** ((self.gamma-1)/self.gamma)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2 / ((self.P_1/self.P_2) ** ((self.gamma-1)/self.gamma))
            return self.T_1

        elif self.P_1 == None:
            self.P_1 = self.P_2 * (self.T_1/self.T_2) ** ((self.gamma-1)/self.gamma)
            return self.P_1

        elif self.P_2 == None:
            self.P_2 = self.P_1 / ((self.T_1/self.T_2) ** ((self.gamma-1)/self.gamma))
            returnself.P_2

        elif self.gamma == None:
            self.gamma = np.log(self.P_2/self.P_1)/(np.log(self.P_2/self.P_1)-np.log(self.T_2/self.T_1))
            return self.gamma


    def gas_law_temp_vol(self):
        if self.T_2 == None:
            self.T_2 = self.T_1 * (self.v_1/self.v_2) ** (self.gamma-1)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2 / ((self.v_1/self.v_2) ** (self.gamma-1))
            return self.T_1

        elif self.P_1 == None:
            self.v_1 = self.v_2 * (self.T_1/self.T_2) ** (self.gamma-1)
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1 / ((self.T_1/self.T_2) ** (self.gamma-1))
            return self.v_2

        elif self.gamma == None:
            self.gamma = (np.log(self.T_2/self.T_1) + np.log(self.v_1/self.v_2)/np.log(self.v_1/self.v_2))
            return self.gamma


    def gamma_cp_cv(self):
        if self.gamma == None:
            self.gamma = self.c_p/self.c_v
            return self.gamma

        elif self.c_p == None:
            self.c_p = self.gamma*self.c_v
            return self.c_p

        elif self.c_v == None:
            self.c_v = self.c_p/self.gamma
            return self.c_v

    def work_energy(self):
        if self.w_12 == None:
            self.w_12 = - ( self.u_2 - self.u_1)
            return self.w_12

        elif self.u_2 == None:
            self.u_2 = self.u_1 - self.w_12
            return self.w_12

        elif self.u_1 == None:
            self.u_1 = self.w_12 - self.u_2
            return self.u_2

    def work_temp(self):
        if self.w_12 == None:
            self.w_12 = -self.c_v*(self.T_2-self.T_1)
            return self.w_12

        elif self.c_v == None:
            self.c_v = -self.w_12/(self.T_2-self.T_1)
            return self.c_v

        elif self.T_2 == None:
            self.T_2 = -self.w_12/self.c_v + self.T_1
            return self.T_2
        
        elif self.T_1 == None:
            self.T_1 = self.T_2 + self.w_12/self.c_v
            return self.T_1
    
    def energy_temp(self):
        if self.u_2 == None:
            self.u_2 = self.u_1 + self.c_v*(self.T_2-self.T_1)
            return self.u_2

        elif self.u_1 == None:
            self.u_1 = self.u_2 - self.c_v*(self.T_2 - self.T_1)
            return self.u_1

        elif self.c_v == None:
            self.c_v = (self.u_2-self.u_1)/(self.T_2-self.T_1)
            return self.c_v

        elif self.T_2 == None:
            self.T_2 = (self.u_2-self.u_1)/self.c_v + self.T_1
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = (self.u_2-self.u_1)/self.c_v - self.T_2
            return self.T_1


    def equation_finder(self):
        called = None
        for i in range(0,8):
            checka = [self.s_1, self.s_2]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Adiabatic.reversible(self)
                called = 1

            checkb = [self.P_1, self.P_2, self.v_1, self.v_2, self.gamma]
            null = 0
            for a in checkb:
                if a == None:
                    null += 1
                    

            if null == 1:
                Adiabatic.gas_law_pres_vol(self)
                called =1

            checkc = [self.T_1, self.T_2, self.P_1, self.P_2, self.gamma]
            null = 0
            for a in checkc:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.gas_law_temp_pres(self)
                called =1


            checkd = [self.T_1, self.T_2, self.v_1, self.v_2, self.gamma]
            null = 0
            for a in checkd:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.gas_law_temp_vol(self)
                called =1

            checke = [self.gamma, self.c_v, self.c_p]
            null = 0
            for a in checkd:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.gamma_cp_cv(self)
                called =1

            checke = [self.w_12, self.u_1, self.u_2]
            null = 0
            for a in checkd:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.work_energy(self)
                called =1

            checkf = [self.w_12, self.T_1, self.T_2, self.c_v]
            null = 0
            for a in checkd:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.work_temp(self)
                called =1

            checkg = [self.u_1, self.u_2, self.T_2, self.c_v, self.T_1]
            null = 0
            for a in checkd:
                if a == None:
                    null += 1

            if null == 1:
                Adiabatic.energy_temp(self)
                called =1



                
        if called == None:
            return [[self.T_1,self.T_2,self.v_1,self.v_2,self.P_1,self.P_2,self.w_12,selfc_v,self.c_p,self.R,self.s_1,self.s_2,self.gamma],['Unable to compute values with variables']]
        else:
            return [[self.T_1,self.T_2,self.v_1,self.v_2,self.P_1,self.P_2,self.w_12,selfc_v,self.c_p,self.R,self.s_1,self.s_2,self.gamma],'N/A']


class Polytropic():
    def __init__(P_1,P_2,v_1,v_2,T_1,T_2,w_12,q_12,c_v,R,n):
        self.P_1 = P_1
        self.P_2 = P_2
        self.v_1 = v_1
        self.v_2 = v_2
        self.T_1 = T_1
        self.T_2 = T_2
        self.w_12 = w_12
        self.q_12 = q_12
        self.c_v = c_v
        self.R = R
        self.n = n

    def gas_law_pres_vol(self):
        if self.P_2 == None:
            self.P_2 = self.P_1 * (self.v_1/self.v_2) ** self.n
            return self.P_2

        elif self.P_1 == None:
            self.P_1 = self.P_2 / ((self.v_1/self.v_2) ** self.n)
            return self.P_1

        elif self.v_1 == None:
            self.v_1 = self.v_2 * (self.P_1/self.P_2) ** self.n
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1 / ((self.P_1/self.P_2) ** ((self.n-1)/self.n))
            self.v_2

        elif self.n == None:
            self.n = np.log(self.P_2/self.P_1)/(np.log(self.v_1/self.v_2))
            return self.n

    def gas_law_temp_pres(self):
        if self.T_2 == None:
            self.T_2 = self.T_1 * (self.P_1/self.P_2) ** ((self.n-1)/self.n)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2 / ((self.P_1/self.P_2) ** ((self.n-1)/self.n))
            return self.T_1

        elif self.P_1 == None:
            self.P_1 = self.P_2 * (self.T_1/self.T_2) ** ((self.n-1)/self.n)
            return self.P_1

        elif self.P_2 == None:
            self.P_2 = self.P_1 / ((self.T_1/self.T_2) ** ((self.n-1)/self.n))
            returnself.P_2

        elif self.n == None:
            self.gamma = np.log(self.P_2/self.P_1)/(np.log(self.P_2/self.P_1)-np.log(self.T_2/self.T_1))
            return self.n


    def gas_law_temp_vol(self):
        if self.T_2 == None:
            self.T_2 = self.T_1 * (self.v_1/self.v_2) ** (self.n-1)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2 / ((self.v_1/self.v_2) ** (self.n-1))
            return self.T_1

        elif self.P_1 == None:
            self.v_1 = self.v_2 * (self.T_1/self.T_2) ** (self.n-1)
            return self.v_1

        elif self.v_2 == None:
            self.v_2 = self.v_1 / ((self.T_1/self.T_2) ** (self.n-1))
            return self.v_2

        elif self.n == None:
            self.n = (np.log(self.T_2/self.T_1) + np.log(self.v_1/self.v_2)/np.log(self.v_1/self.v_2))
            return self.n    

    def work_pres_vol(self):
        if self.w_12 == None:
            self.w_12 = (self.P_2*self.v_2 - self.P_1*self.v_1)/(1-self.n)
            return self.w_12

        elif self.P_2 == None:
            self.P_2 = (self.w_12*(1-self.n) + self.P_1*self.v_1)/self.v_2
            return self.P_2

        elif self.P_1 == None:
            self.P_1 = (self.w_12*(1-self.n) - self.P_2*self.v_2)/self.v_1
            return self.P_1

        elif self.v_2 == None:
            self.v_2 = (self.w_12*(1-self.n) + self.P_1*self.v_1)/self.P_2
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = (self.w_12*(1-self.n) - self.P_2*self.v_2)/self.P_1
            return self.v_2

        elif self.n == None:
            self.n = 1 - (self.P_2*self.v_2 - self.P_1*self.v_1)/(self.w_12)
            return self.w_12

    def work_temp_R(self):
        if self.w_12 == None:
            self.w_12 = self.R*(self.T_2 - self.T_1)/(1-self.n)
            return self.w_12

        elif self.R == None:
            self.R = self.w_12*(1-self.n)/(self.T_2-self.T_1)
            return self.R

        elif self.T_2 == None:
            self.T_2 = self.w_12*(1-self.n)/self.R + self.T_1
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = - self.w_12*(1-self.n)/self.R + self.T_2
            return self.T_1

        elif self.n == None:
            self.n = 1 - self.R/self.w_12*(self.T_2 - self.T_1)
            return self.n

    def pres_vol_temp(self):
        if self.P_2 == None:
            self.P_2 = (self.R*(self.T_2-self.T_1) + self.P_1+self.v_1)/self.v_2
            return self.P_2

        elif self.P_1 == None:
            self.P_1 = (-self.R*(self.T_2-self.T_1) + self.P_2+self.v_2)/self.v_1
            return self.P_1

        elif self.v_2 == None:
            self.v_2 = (self.R*(self.T_2-self.T_1) + self.P_1+self.v_1)/self.P_2
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = (-self.R*(self.T_2-self.T_1) + self.P_2+self.v_2)/self.P_1
            return self.v_1

        elif R == None:
            self.R = (self.P_2*self.v_2 - self.P_1*self.v_1)/(self.T_2-self.T_1)
            return self.R

        elif T_2 == None:
            self.T_2 = (self.P_2*self.v_2 - self.P_1*self.v_1)/(self.R) + self.T_1
            return self.T_2

        elif T_1 == None:
            self.T_1 = -(self.P_2*self.v_2 - self.P_1*self.v_1)/(self.R) + self.T_2
            return self.T_1

    def heat_transferred(self):
        if self.q_12 == None:
            self.q_12 = (self.c_v + self.R/(1-self.n))(self.T_2-self.T_1)
            return self.q_12

        elif self.c_v == None:
            self.c_v = self.q_12/(self.T_2-self.T_1)-self.R/(1-self.n)
            return self.c_v

        elif self.R == None:
            self.R = (1-self.n)*(self.q_12/(self.T_2-self.T_1-self.c_v))
            return self.R

        elif self.n == None:
            self.n = 1-1/((self.q_12/(self.T_2-self.T_1)-self.c_v)/self.R)
            return self.n

        elif self.T_2 == None:
            self.T_2 = self.T_1 + self.q_12/(self.c_v+self.R/(1-self.n))
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2 - self.q_12/(self.c_v+self.R/(1-self.n))
            return self.T_1


    def entropy_temp_vol(self):
        if self.s_2 == None:
            self.s_2 = self.s_1 + self.c_v*np.log(self.T_2/self.T_1) + self.R*np.log(self.v_2/self.v_1)
            return self.s_2
        
        elif self.s_1 == None:
            self.s_1 = self.s_2 - self.c_v*np.log(self.T_2/self.T_1) + self.R*np.log(self.v_2/self.v_1)
            return self.s_1

        elif self.c_v == None:
            self.c_v = (self.s_2 - self.s_1 - self.R*np.log(self.v_2/self.v_1))/np.log(self.T_2/self.T_1)
            return self.c_v
        
        elif self.T_2 == None:
            self.T_2 = self.T_1*np.exp((self.s_2 - self.s_1 - self.R*np.log(self.v_2/self.v_1))/self.c_v)
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2/np.exp((self.s_2 - self.s_1 - self.R*np.log(self.v_2/self.v_1))/self.c_v)
            return self.T_1
        
        elif self.v_2 == None:
            self.v_2 = self.v_1*np.exp((self.s_2 - self.s_1 - self.c_v*np.log(self.T_2/self.T_1))/self.R)
            return self.v_2

        elif self.v_1 == None:
            self.v_1 = self.v_2/np.exp((self.s_2 - self.s_1 - self.c_v*np.log(self.T_2/self.T_1))/self.R)
            return self.v_1

    
    def entropy_temp(self):
        if self.s_2 == None:
            self.s_2 = self.s_1 + (self.c_v + self.R(1-self.n))*np.log(self.T_2/self.T_1)
            return self.s_2

        elif self.s_1 == None:
            self.s_1 = self.s_2 - (self.c_v + self.R(1-self.n))*np.log(self.T_2/self.T_1)
            return self.s_1

        elif self.c_v == None:
            self.c_v = (self.s_2 - self.s_1)/np.log(self.T_2/self.T_1) - self.R/(1 - self.n)
            return self.c_v

        elif self.R == None:
            self.R = (1-self.n)*((self.s_2 - self.s_1)/np.log(self.T_2/self.T_1) - self.c_v)
            return self.R

        elif self.n == None:
            self.n = 1 -1/(((self.s_2-self.s_1)/np.log(self.T_2/self.T_1)-self.c_v)/self.R)
            return self.n

        elif self.T_2 == None:
            self.T_2 = self.T_1*np.exp((self.s_2 - self.s_1)/(self.c_v + self.R/(1-self.n)))
            return self.T_2

        elif self.T_1 == None:
            self.T_1 = self.T_2/np.exp((self.s_2 - self.s_1)/(self.c_v + self.R/(1-self.n)))
            return self.T_1



    def equation_finder(self):
        called = None
        for i in range(0,9):
            checka = [self.P_1, self.P_2, self.v_1, self.v_2, self.n]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.gas_law_pres_vol(self)
                called = 1

            checkb = [self.P_1, self.P_2, self.T_1, self.T_2, self.n]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.gas_law_temp_pres(self)
                called = 1

            checkc = [self.T_1, self.T_2, self.v_1, self.v_2, self.n]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.gas_law_temp_vol(self)
                called = 1


            checkd = [self.w_12, self.P_1, self.P_2, self.v_1, self.v_2, self.n]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.work_pres_vol(self)
                called = 1

            checke = [self.w_12, self.R, self.T_2, self.T_1, self.n]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.work_temp_R(self)
                called = 1

            checke = [self.w_12, self.P_1, self.P_2, self.v_2, self.v_1, self.T_2, self.T_1]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.pres_vol_temp(self)
                called = 1

            checkf = [self.q_12, self.R, self.n, self.T_1, self.T_2, self.c_v]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.heat_transferred(self)
                called = 1

            checkg = [self.s_1, self.s_2, self.T_2, self.c_v, self.T_1, self.v_2, self.v_1]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.entropy_temp_vol(self)
                called = 1


            checkh = [self.s_1, self.s_2, self.T_2, self.c_v, self.T_1, self.n, self.R]
            null = 0
            for a in checka:
                if a == None:
                    null += 1
            if null == 1:
                Polytropic.entropy_temp(self)
                called = 1




                
        if called == None:
            return [[self.P_1,self.P_2,self.v_1,self.v_2,self.T_1,self.T_2,self.w_12,self.q_12,self.c_v,self.R,self.n],['Unable to compute values with variables']]
        else:
            return [[self.P_1,self.P_2,self.v_1,self.v_2,self.T_1,self.T_2,self.w_12,self.q_12,self.c_v,self.R,self.n],'N/A']

        
if __name__ == "__main__":

    constVol = Constant_Volume(123,23,43,None,None,None,None,None)
    a = constVol.equation_finder()
    print(a)
    constPres = Constant_Pressure(50,None,300,200,5,None,None,4000,None,400,300,30)
    b = constPres.equation_finder()
    print(b)
    #T_1,T_2,v_1,v_2,q_12,w_12,c_v,c_p,R,s_1,s_2,p
    constTemp = Constant_Temperature(1,2,3,4,5,6,7,8,9,10,11,12)
    b = constTemp.equation_finder()
    print(b)

    #class Constant_Temperature():
    #def __init__(self,T,v_1,v_2,P_1,P_2,q_12,w_12,c_v,c_p,R,s_1,s_2):

    