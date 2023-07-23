import numpy as np
import openseespy.opensees as ops
from scipy import optimize
from scipy.stats import norm
from pathlib import Path
import opsvis as opsv

# To Do
# 1.- improve the set_support methos it should be a dict or a list of list 

class Model():
    def __init__(self, node_path, conectivity_path, element_id_path, members_prop_path, members_name_path, mesh_cells_path, surfaces_path) -> None:
        self.nodes = np.genfromtxt(node_path, invalid_raise=False)
        self.conect = np.loadtxt(conectivity_path, delimiter=',')
        self.idele = np.loadtxt(element_id_path, delimiter=',')
        self.members_prop = np.genfromtxt(members_prop_path, invalid_raise=False)
        self.names = np.genfromtxt(members_name_path, dtype='str', delimiter = '\t')
        self.mesh_cells = np.genfromtxt(mesh_cells_path, invalid_raise=False)
        self.surfaces = np.genfromtxt(surfaces_path, invalid_raise=False)
        self.nodes_w_support = None
        self.default_restrains = None

    @property
    def steel_properties(self):
        G = 0.25*10**3
        gama = 78.5*10**-6
        E = 200*10**3 
        return G,gama,E
    
    @property
    def concrete_properties():
        G_c = 11.125 # N/mm2
        gama_conc = 25*10**-6
        E_conc = 26 # N/mm2
        
        return G_c, gama_conc, E_conc
        

    def set_supports(self,nodes_w_support,default_restrains = [1,1,1,1,1,1]):
        self.nodes_w_support = nodes_w_support
        self.default_restrains = default_restrains

    def create_nodes_mass(self):  
        # Create a numpy array to alocate mass 
        masa = np.zeros((len(self.nodes),2))
        # Create opensees nodes 
        for i in range(len(self.nodes)):
            ops.node(i+1,self.nodes[i,0],self.nodes[i,1],self.nodes[i,2])
            masa[i, 0] = i+1
        return masa
    
    def internal_Forces(self, i, E, G_mod, N):
        #This function creates: 
        #the recorders for each fucntions
        #the elastic section 
        #section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
        #the integrations points with Lobatos
        #beamIntegration('Lobatto', tag, secTag, N)
    
        A = self.members_prop[int(self.idele[i])-1,3]
        Iz = self.members_prop[int(self.idele[i])-1,1]
        Iy = self.members_prop[int(self.idele[i])-1,2]
        Jxx = self.members_prop[int(self.idele[i])-1,0]
        
        ops.section('Elastic',i+1, E, A, Iz, Iy, G_mod, Jxx)
        ops.beamIntegration('Lobatto', i+1, i+1, N)  

    
    def create_mass_elements(self, masa, z, gama, g, E, G, N, verbose=False):
        
        for i in range(len(self.conect[:,1])):
            
            # Node i&j Coordinate
            XYZI = ops.nodeCoord(int(self.conect[i,0]))
            XYZJ = ops.nodeCoord(int(self.conect[i,1]))
            # xaxis is parales to the element axis
            xaxis = np.subtract(XYZJ,XYZI)  
            vecxz = np.cross(xaxis,z)
            # Length of elements
            L =np.linalg.norm(xaxis)
            # Mass for each sections
            m =self.members_prop[int(self.idele[i])-1,3]*L*gama/2/g
        
        
            if np.linalg.norm(vecxz) == 0:
        
                ops.geomTransf('Linear', i,0,-1,0)
                self.internal_Forces(i, E, G, N) 
                con = [int(self.conect[i,0]),int(self.conect[i,1])]       
                ops.element('forceBeamColumn',i+1,*con,i,i+1)
                index = np.where(masa[:,0] == self.conect[i,0])[0][0]
                masa[index,1] += m
                index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                masa[index,1] += m
            else:
                
                if XYZI[2]-XYZI[2] == 0:
        
                    ops.geomTransf('Linear', i,*vecxz)
                    self.internal_Forces(i,E,G,N)
                    con = [int(self.conect[i,0]),int(self.conect[i,1])]       
                    ops.element('forceBeamColumn',i+1,*con,i,i+1)
                    index = np.where(masa[:,0] == self.conect[i,0])[0][0]
                    masa[index,1] += m
                    index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                    masa[index,1] += m
                else:
                    #posible bug ?
                    #self.internal_Forces(i,E,G,N)
                    ops.geomTransf('Linear', i,*vecxz)
                    self.internal_Forces(i,E,G,N)
                    con = [int(self.conect[i,0]),int(self.conect[i,1])]       
                    ops.element('Truss',i+1,*con,i,i+1)
                    index = np.where(masa[:,0] == self.conect[i,0])[0][0]
                    masa[index,1] += m
                    index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                    masa[index,1] += m

        if verbose:
            print('Section  ',self.names[int(self.idele[i])-1],'  assigned to element No.  ',i+1)

    def fix_model(self):   
        for i in range(len(self.nodes_w_support)):
            ops.fix(int(self.nodes_w_support[i]),*self.default_restrains)

    def assign_mass(self,masa):
        for i in range(len(masa)):
             
             ops.mass(masa[i, 0], masa[i, 1], masa[i, 1], masa[i, 1], 0, 0, 0)

    def Modal_analysis(Nm,sw):
    
        """
        Nm: Number of Modes that the code is going to return 
        sw: Switch value for returning coordinates 
        """
        Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150]
        Const1 = [17,70,138,133,128,97,37,7,12]
        freq = ops.eigen(Nm)
        Nodes_AC= np.concatenate((Const1,Const2))
        
        U_ops = np.zeros((36,Nm))
        for i in range(0,Nm):
            for j in range(0,len(Nodes_AC)*2):
                if j<18:
                    U_ops[j,i]= ops.nodeEigenvector(int(Nodes_AC[j]),i+1)[0]
                else:
                    U_ops[j,i]= ops.nodeEigenvector(int(Nodes_AC[j-18]),i+1)[1]
        
        
        for i in range(0,Nm):
            U_ops[:,i] =  U_ops[:,i]/np.linalg.norm(U_ops[:,i])
        
        # print('Mode:','Frequency [Hz]')
        # for i in range(Nm):
        #     print(i+1,freq[i]**0.5/(2*np.pi))
        
        if sw == 1:
            CORR = np.zeros((len(Nodes_AC),4))
            for i  in range(len(Nodes_AC)):
                CORR[i,0] = Nodes_AC[i]
                CORR[i,1:5] = ops.nodeCoord(int(Nodes_AC[i]))
                print(CORR[i,:])
                return freq,U_ops,CORR
        else:
            CORR =[]
            return freq,U_ops
                        
    def create_model(self,verbose,Nm):
        ops.wipe()
        # Creates Opensees Model
        ops.model('basic','-ndm',3,'-ndf',6)
        # Z global axis reference
        z = [0,0,1]
        # integration points 
        N = 6
        # generate nodes and 
        G,gama,E = self.steel_properties
        # define gravity 
        g = 9800
        # Assign Mases 
        masa = self.create_nodes_mass()

        self.create_mass_elements(masa, z, gama, g, E, G, N, verbose)
        # assign self weight masses 
        self.assign_mass(masa)
        self.fix_model()
        freq = ops.eigen(Nm)

        if verbose:
            print('Current Model Frequencies')
            for i in range(len(freq)):
                print(freq[i]**0.5/(2*np.pi))


        

node_path = 'Nodes.txt' 
conectivity_path = 'Connc.txt'
element_id_path = 'Elementid.txt'
members_prop_path = 'Members.txt'
members_name_path = 'Nmes.txt'

support = [1,4,9,14,29,34,65,67,89,94,122,125,130,135,139,142,147\
            ,152,170,175,207,209,230,235,260,263,268,273]

model = Model(node_path,conectivity_path,element_id_path,members_prop_path,members_name_path)
model.set_supports(support)
model.create_model(verbose = True ,Nm =  10)
opsv.plot_mode_shape(3,1000)
