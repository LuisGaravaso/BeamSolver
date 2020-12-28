# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 14:35:50 2020

@author: otavi
"""

#Imports
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from functools import partial

#Setup
sns.set()
colors = sns.color_palette()

#Classes
class Sing:
    
    def __init__(self, mag, pos, expo):
        
        '''
        The singularity is a function that works like a switch.
        It is written as <x - pos>**expo.
        It evaluates into '0' if (x < pos) or if (expo < 0).
        '''
        
        self.mag = mag #Magnitude
        self.pos = pos #Position of the Singularity
        self.expo = expo #Exponent of the Singularity
        self._is_supp = False #Support Flag used by the Supp() class
        
    def evaluate(self,x):
        
        '''
        Função que retorna o valor numérico da Singularidade
        
        Parâmetros
        ----------
        x : float
            coordenada 'x' que passa pela viga 
    
        Returna
        -------
        valor numérico da singularidade: float
        '''
        
        mag = self.mag
        pos = self.pos
        expo = self.expo
        
        if (x < pos) or (expo < 0):
            return 0
        elif isinstance(expo, float):
            raise ValueError("O expoente deve ser um número inteiro")
            
        else:
            return mag*(x - pos)**expo
    
    def integrate(self):
        '''
        
        Função que faz a integral da singularidade
        
        Returns
        -------
        Sing()
            Retorna uma nova singularidade integrada
        '''
        
        mag = self.mag
        pos = self.pos
        expo = self.expo
        
        n = 1
        if expo >= 1:
            n = expo + 1
        
        mag_ = mag/n
        expo_ = expo + 1
        S = Sing(mag_, pos, expo_)
        S._is_supp = self._is_supp
        return S
    
    def __add__(self, other):
        
        sings = [self, other]
        return Carregamento(sings)
    
    def __repr__(self):
        
        mag = self.mag
        pos = self.pos
        expo = self.expo
        
        if self._is_supp:
            if mag != 1:
                mag = f"{mag:.4f}"+"⋅Ry"
            else:
                mag = "Ry"
            
            return f"{mag} ⋅ <x - {pos}>^({expo})"
            
        if (mag == 0):
            strmaker = ""
        elif (pos == 0) and (expo == 0):
            strmaker = f"{mag:.4f}"
        elif pos > 0:
            strmaker = f"{mag:.4f} ⋅ <x - {pos}>^({expo})"
        else:
            strmaker = f"{mag:.4f} ⋅ x^({expo})"
        
        return strmaker

def Supp(pos):
    
    '''
    Support Singularity.
    This class flags Sing()._is_supp as True.
    At first the Reaction on the Support is unknown, so it must be
    evaluated by the solver.
    
    Parameters
    ----------
    pos : float
        position of the Support Singularity

    Returna
    -------
    Supp: Sing()
        The Support singularity
    '''
    
    mag = 1
    pos = pos
    expo = -1
    
    S = Sing(mag, pos, expo)
    S._is_supp = True
    
    return S
        
class Carregamento:
    
    def __init__(self, Sings):
        '''
        Contrutor:
            Recebe Lista de Singularidades
        '''
        Sings.sort(key = lambda x: x.pos)
        
        PosLst = [sing.pos for sing in Sings]
        Positions = sorted(set([0,*PosLst]))
        Mapper = {pos: chr(97 + i).upper() for i, pos in enumerate(Positions)}
        
        self.Mapper = Mapper
        self.Sings = Sings
        self.__PosLst = PosLst
        
    def evaluate(self, x, force = False):
        
        Apoio = []
        Cargas = []
        for sing in self.Sings:
            value = sing.evaluate(x)
            if sing._is_supp:
                Apoio.append(value)
            else:
                Cargas.append(value)
        
        if force:
            return sum(Cargas) + sum(Apoio)
    
        return Cargas, Apoio
    
    def integrate(self):
        
        Sings_ = []
        for sing in self.Sings:
            Sings_.append(sing.integrate())
        
        return Carregamento(Sings_)
        
    def __add__(self, other):
        
        sings = self.Sings
        if isinstance(other, Sing):
            sings.append(other)
        
        elif isinstance(other, Carregamento):
            sings.extend(other.Sings)
    
        else:
            raise ValueError("Os objetos devem ser Sing() ou Carregamento()")
            
        return Carregamento(sings)
    
    def __repr__(self,):
        strmaker = ""
        for pos, sing in zip(self.__PosLst, self.Sings):
            label = self.Mapper[pos]
            text = sing.__repr__().replace("Ry",f"Ry{label}")
            
            if text == "":
                strmaker += ""
            else:
                strmaker += text + " + "
        
        if strmaker == "":
            return "0"
        else:
            return strmaker[:-3]
    

class Viga:
    
    def __init__(self, Q, L, E = None, Izz = None):
        self.Q = Q #Equação do Carregamento
        self.L = L#Comprimento Total da Viga
        self.E = E #Módulo de Elasticidade do Material
        self.Izz = Izz #Momento de Área da Seção Transversal
        self.Conds = defaultdict(list) #Condições de Contorno e Restrição
        self._solved = False #Flag que verifica se o problema foi resolvido
        
        if Q == 0: #Caso de Carregamento Nulo
            self.Q = Carregamento([Sing(0,0,1)])
        
        elif isinstance(Q, Sing):
            self.Q = Carregamento([Q])
        
        self.Mapper = self.Q.Mapper
        letter = list(self.Mapper.values())[-1]
        letter = chr(ord(letter) + 1).upper()
        self.Mapper[self.L] = letter
            
        self._NCond = 0
        if (E is None) or (Izz is None):
            self.NCond = 2
        else:
            self.NCond = 4
        
        for sing in self.Q.Sings:
            if (sing._is_supp) and (self.NCond == 4):
                pos = sing.pos
                self.Conds[pos].append(("Flecha", 0))    
        
    def add_forca(self, mag, pos):
        add = self.__assertConds(flag = 0)
        pos = self.__assertpos(pos)
        
        if add:
            self.Conds[pos].append(("V",mag))
        
    def add_momento(self, mag, pos):
        add = self.__assertConds(flag = 0)
        pos = self.__assertpos(pos)
        
        if add:
            self.Conds[pos].append(("M",mag))
    
    def add_engaste(self, pos):
        add = self.__assertConds(flag = 1)
        pos = self.__assertpos(pos)
        
        if add:
            self.Conds[pos].append(("Theta", 0))
            self.Conds[pos].append(("Flecha", 0))
        
    def add_supp(self, pos):
        add = self.__assertConds(flag = 2)
        pos = self.__assertpos(pos)
        
        if add:
            self.Conds[pos].append(("Flecha", 0))
    
    def solve(self):
        
        #Constantes de Integração 
        C1 = Sing(1, 0, 0)
        C2 = Sing(1, 0, -1)
        Cs = C1 + C2
        
        if self.NCond > 2:
            C3 = Sing(1, 0, -2)
            C4 = Sing(1, 0, -3)
            Cs += C3 + C4
        
        #Integração das Constantes
        Cs_1l = Cs.integrate()
       
        if self.NCond > 2:
            Cs_2l = Cs_1l.integrate()
            Cs_3l = Cs_2l.integrate()
        
        #Matrizes do Sistema
        A = []
        b = []
        
        #Equações de Esforços
        Q = self.Q
        V = Q.integrate()
        M = V.integrate()
        EITeta = M.integrate()
        EIFlecha = EITeta.integrate() 
        
        #Valores Numéricos
        for pos, conds in self.Conds.items():
            for tipo, value in conds:
                if tipo == "V":
                    cargas, apoios = V.evaluate(pos)
                    ctes, _ = Cs.evaluate(pos)
                elif tipo == "M":
                    cargas, apoios = M.evaluate(pos)
                    ctes, _ = Cs_1l.evaluate(pos)
                elif tipo == "Theta":
                    cargas, apoios = EITeta.evaluate(pos)
                    ctes, _ = Cs_2l.evaluate(pos)
                else:
                    cargas, apoios = EIFlecha.evaluate(pos)
                    ctes, _ = Cs_3l.evaluate(pos)
                
                ctes.extend(apoios)
                A.append(ctes)
                b.append(value - sum(cargas))
        
        X = np.linalg.solve(A,b)
        X = X.round(8)
        
        if self.NCond == 2:
            if len(X) > 2:
                C1, C2, *Apoios = X
            else:
                C1, C2 = X
                Apoios = []
        
        else:
            if len(X) > 4:
                C1, C2, C3, C4, *Apoios = X
            else:
                C1, C2, C3, C4 = X
                Apoios = []
        
        #Atualizar Constantes de Integração 
        C1 = Sing(C1, 0, -1)
        C2 = Sing(C2, 0, -2)
        Cs = C1 + C2
        
        if self.NCond > 2:
            C3 = Sing(C3, 0, -3)
            C4 = Sing(C4, 0, -4)
            Cs += C3 + C4
        
        #Atualizar Apoios
        i = 0
        for ind, sing in enumerate(self.Q.Sings):
            if sing._is_supp:
                pos = sing.pos
                mag = Apoios[i]
                self.Q.Sings[ind] = Sing(mag, pos, -1)
                i += 1
                
        self.Q = self.Q + Cs
        
        #Equações de Carregamento
        self.Vy = self.Q.integrate()
        self.Mz = self.Vy.integrate()
        
        if self.NCond > 2:
            self.EITheta = self.Mz.integrate()
            self.EIFlecha = self.EITheta.integrate()
        
        self._solved = True
        
    def plot(self,):
        
        if not self._solved:
            raise ValueError("O método Viga.solve() ainda não foi chamado")
        
        Eqs = [self.Vy, self.Mz]
        Titles = ["$V_y(x)$", "$M_z(x)$"]
        if self.NCond == 2:
            fig, axs = plt.subplots(2,1, sharex = True, sharey = True, figsize = (8,16))
        else:
            fig, axs = plt.subplots(2,2, sharex = True, sharey = 'row', figsize = (12,12))
            Eqs.append(self.EITheta)
            Eqs.append(self.EIFlecha)
            Titles.append(r"$EI \theta_z (x)$")
            Titles.append(r"$EI \nu (x)$")

        i = 0
        X = np.linspace(0, self.L, num = 1000)
        Y0 = np.zeros(X.shape)
        
        for ax, Eq, title in zip(axs.flatten(), Eqs, Titles):
            
            func = partial(Eq.evaluate, force = True)
            func = np.vectorize(func)
            Y = func(X)
            
            ax.hlines(0, 0, self.L, color = "black")
            ax.plot(X, Y, label = title, color = colors[i])
            ax.fill_between(X, Y0, Y, color = colors[i], alpha = 0.5)
            ax.set_xlabel("$x\ [m]$", fontsize = 14)
            ax.set_ylabel(title, fontsize = 14)
            ax.tick_params(labelbottom=True)
            ax.tick_params(labelleft=True)
            ax.legend()
            
            ymin, ymax = ax.get_ylim()
            for pos, value in self.Mapper.items():
                ax.vlines(pos, ymin, ymax, linestyle = "dashed", alpha = 0.7)
                
                if pos < self.L:
                    ax.annotate(value,(pos + (self.L/100),0.95*ymin))
                else:
                    ax.annotate(value,(pos - 4*(self.L/100),0.95*ymin))
                    
            ax.set_ylim(ymin, ymax)
            
            i += 1
            
        ax.set_xlim([-0.05,self.L + 0.05])
        return fig
    
    def __assertpos(self, pos):
        if (pos == "i") or (pos == 0):
            return 0
        else:
            return self.L
    
    def __assertConds(self, flag):
        E = self.E
        Izz = self.Izz
        add = False
        if flag == 0: #Forca ou Momento
            self._NCond += 1
            add = True
            
        elif flag == 1: #Engaste
            if (E is None) or (Izz is None):
                self._NCond += 0
            else:
                self._NCond += 2
                add = True
                
        else: #Apoio
            if (E is None) or (Izz is None):
                self._NCond += 0 
            else:
                self._NCond += 1
                add = True
                
        if self._NCond < self.NCond:
            print(f"{self.NCond - self._NCond} condição/condições faltantes")
            
        elif self._NCond == self.NCond:
            print("Número de Condições Satisfeitas. Chamar Viga.solve()")
        
        else:
            raise ValueError("Foram fornecidas mais condições do que necessário")
        
        return add
    
print("Este Solver requer a utilização da seguinte convenção de sinais:")
conv = plt.imread("Convencao.png")
plt.imshow(conv, aspect = "equal")
plt.axis('off')
print("Fonte: Apostila de Resistência dos Materiais I")
print("https://sites.google.com/view/jlabaki/teaching/resmat-1")
