import numpy as np
import openseespy.opensees as ops
from scipy import optimize
from scipy.stats import norm
from pathlib import Path
import opsvis as opsv

# To Do
# 1.- improve the set_support method it should be a dict or a list of list 

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
        self.global_surface_tag = 100000.0
        self.g = 9800
    @property
    def steel_properties(self):
        G = 0.25*10**3
        gama = 7.85*10**-5
        E = 200*10**3 
        return G,gama,E
    

    def concrete_properties(self):
        G_c = 11125 # N/mm2
        gama_conc = 2.5*10**-5
        E_conc = 26700  # N/mm2
        
        return G_c, gama_conc, E_conc
        

    def get_index_cross_section_prop(self, _member_id):
        # idele tells me what cross section number the membe id has assinged
        idele_location = np.where(_member_id == self.idele[:,0])[0][0]
        _idele = self.idele[idele_location,1]
        rotation = self.idele[idele_location,2]
        # we het the cross section number now we look into cross section props
        # to get the index of the cross section number 
        cross_index = np.where(_idele == self.members_prop[:,0])[0][0]
        return cross_index , rotation
        
        

    def set_supports(self,nodes_w_support,default_restrains = [1,1,1,1,1,1]):
        self.nodes_w_support = nodes_w_support
        self.default_restrains = default_restrains

    def create_nodes_mass(self):  
        # Create a numpy array to alocate mass 
        masa = np.zeros((len(self.nodes),2))
        # Create opensees nodes 
        for i in range(len(self.nodes)):
            #id,X,Y,Z
            ops.node(self.nodes[i,0],self.nodes[i,1],self.nodes[i,2],self.nodes[i,3])
            masa[i, 0] = self.nodes[i,0]
        return masa

    def mass_shell(self, i, j, k, l, gama_conc, t):
        # assuming ops.nodeCoord returns a list of coordinates [x, y, z]
        cords_i = np.array(ops.nodeCoord(int(i)))
        cords_j = np.array(ops.nodeCoord(int(j)))
        cords_k = np.array(ops.nodeCoord(int(k)))
        cords_l = np.array(ops.nodeCoord(int(l)))
        
        # define vectors of the triangles
        vector_ij = cords_j - cords_i
        vector_ik = cords_k - cords_i
        vector_il = cords_l - cords_i
        
        # calculate areas of the two triangles (1/2 * magnitude of the cross product of two sides)
        area1 = np.linalg.norm(np.cross(vector_ij, vector_ik)) / 2.0
        area2 = np.linalg.norm(np.cross(vector_ij, vector_il)) / 2.0
        
        # total area of the quadrilateral
        area = area1 + area2
        
        # calculate mass
        mass = t * area * gama_conc/self.g/4
        return mass
    
    def create_shells(self, masa, verbose=False):
        surface_tag = self.surfaces[0]
        # thickness
        t = self.surfaces[2] 
        G_c, gama_conc, E_conc = self.concrete_properties()
        gtag = self.global_surface_tag
        ops.section('ElasticMembranePlateSection', int(surface_tag+gtag), E_conc, 0.2, t, 0)
        for shell_id, _, i, j, k, l in self.mesh_cells:
            # check gtag jejej
            idd = shell_id+gtag
            ops.element('ShellMITC4',int(idd), int(i), int(j), int(k), int(l), int(surface_tag+gtag))
            # Get Mass
            m= self.mass_shell(i, j, k, l, gama_conc, t)
            
            index = np.where(masa[:,0] == int(i))[0][0]
            masa[index,1] += m
            index = np.where(masa[:,0] == int(j))[0][0]
            masa[index,1] += m
            index = np.where(masa[:,0] == int(k))[0][0]
            masa[index,1] += m
            index = np.where(masa[:,0] == int(l))[0][0]
            masa[index,1] += m
            if verbose:
                pass#print(f"the id :{shell_id}, with the nodes 1 :{i} , 2 : {j} , 3:{k}, 4:{l} @ with mass : {m}")
                
                
        return masa 
        
    
    """ Requires Modification """ 
    def internal_Forces(self, _member_id , cross_index ,rotation, E, G_mod, N):
        #This function creates: 
        #the recorders for each fucntions
        #the elastic section 
        #section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
        #the integrations points with Lobatos
        #beamIntegration('Lobatto', tag, secTag, N)
        #"""""""""""""""""""""""""""""""""""""""""""""#
       
        if rotation != 0:
            A = self.members_prop[cross_index,4]
            Iz = self.members_prop[cross_index,3]
            Iy = self.members_prop[cross_index,2]
            Jxx = self.members_prop[cross_index,1]
            print('this member {self.names[cross_index]} was rotated 90')
        else:
            
            A = self.members_prop[cross_index,4]
            Iz = self.members_prop[cross_index,2]
            Iy = self.members_prop[cross_index,3]
            Jxx = self.members_prop[cross_index,1]            
        
        ops.section('Elastic',_member_id, E, A, Iz, Iy, G_mod, Jxx)
        ops.beamIntegration('Lobatto', _member_id, _member_id, N)  
        
    
    def create_mass_elements(self, masa, z, gama, g, E, G, N, verbose=False):
        
        for i in range(len(self.conect[:,1])):
            # get member id 
            members_id = int(self.conect[i,0])
            # get cross secction associate to member id 
            cross_index,rotation = self.get_index_cross_section_prop(members_id);
                  
            # Node i&j Coordinate
            XYZI = ops.nodeCoord(int(self.conect[i,1]))
            XYZJ = ops.nodeCoord(int(self.conect[i,2]))
            # xaxis is parales to the element axis
            xaxis = np.subtract(XYZJ,XYZI)  
            vecxz = np.cross(xaxis,z)
            # Length of elements
            L =np.linalg.norm(xaxis)
            # Mass for each sections
            m =self.members_prop[cross_index,4]*L*gama/2/g
        
        
            if np.linalg.norm(vecxz) == 0:
        
                ops.geomTransf('Linear',members_id ,0,-1,0)
                self.internal_Forces(members_id , cross_index,rotation , E, G, N) 
                con = [int(self.conect[i,1]),int(self.conect[i,2])]   
                #element('forceBeamColumn', eleTag, *eleNodes, transfTag,
                #integrationTag, '-iter', maxIter=10, tol=1e-12, '-mass', mass=0.0)
                ops.element('forceBeamColumn',members_id,*con,members_id,members_id)
                index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                masa[index,1] += m
                index = np.where(masa[:,0] == self.conect[i,2])[0][0]
                masa[index,1] += m
            else:
                
                if XYZI[2]-XYZI[2] == 0:
        
                    ops.geomTransf('Linear', members_id ,*vecxz)
                    self.internal_Forces(members_id ,cross_index,rotation, E,G,N)
                    con = [int(self.conect[i,1]),int(self.conect[i,2])]       
                    ops.element('forceBeamColumn',members_id ,*con,members_id ,members_id )
                    index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                    masa[index,1] += m
                    index = np.where(masa[:,0] == self.conect[i,2])[0][0]
                    masa[index,1] += m
                else:
                    #posible bug ?
                    #self.internal_Forces(i,E,G,N)
                    ops.geomTransf('Linear', members_id,*vecxz)
                    self.internal_Forces(members_id,cross_index,rotation ,E,G,N)
                    con = [int(self.conect[i,1]),int(self.conect[i,2])]       
                    ops.element('Truss',members_id,*con,members_id,members_id)
                    index = np.where(masa[:,0] == self.conect[i,1])[0][0]
                    masa[index,1] += m
                    index = np.where(masa[:,0] == self.conect[i,2])[0][0]
                    masa[index,1] += m

        if verbose:
            #print('Section  ',self.names[int(self.idele[i])-1],'  assigned to element No.  ',i+1)
            print('To complex')

    def fix_model(self):   
        for i in range(len(self.nodes_w_support)):
            ops.fix(int(self.nodes_w_support[i]),*self.default_restrains)

    def assign_mass(self,masa):
        for i in range(len(masa)):
             ops.mass(masa[i, 0], masa[i, 1], masa[i, 1], masa[i, 1], 0, 0, 0)

    def Modal_analysis(self, Nm):
        freq = ops.eigen(Nm)
        return freq
    
                        
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
        g = self.g
        # Assign Mases 
        masa = self.create_nodes_mass()

        self.create_mass_elements(masa, z, gama, g, E, G, N, verbose)
        #crate mass 
        if self.mesh_cells.size ==0:
            
            self.assign_mass(masa)
            print('masas asigned')
        else:
            masa = self.create_shells( masa, verbose)
            self.assign_mass(masa)         

        # assign self weight masses 
        
        self.fix_model()
        freq = self.Modal_analysis(Nm)

        if verbose:
            print('Current Model Frequencies')
            for i in range(len(freq)):
                print(freq[i]**0.5/(2*np.pi))
        return ops , masa 


        

# node_path = 'Nodes.txt' 
# conectivity_path = 'Connc.txt'
# element_id_path = 'Elementid.txt'
# members_prop_path = 'Members.txt'
# members_name_path = 'Nmes.txt'

# support = [1,4,9,14,29,34,65,67,89,94,122,125,130,135,139,142,147\
#             ,152,170,175,207,209,230,235,260,263,268,273]

# model = Model(node_path,conectivity_path,element_id_path,members_prop_path,members_name_path)
# model.set_supports(support)
# model.create_model(verbose = True ,Nm =  10)
# opsv.plot_mode_shape(3,1000)
