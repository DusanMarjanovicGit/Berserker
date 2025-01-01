import inca 
import inca.component.geometry.closegeometry as cg 
from inca.component.geometry.tube.cylinder import Cylinder 
import sys 
sys.path.append("/home/dm266236/AP3/StaticCase/Coeur/Minaret/3D/Robustesse/General/ConceptionPackage")
import MeshManager
import SectionSetter
import numpy  as np
import ShapeSlicer
from math import pi 

onMaterial="MATERIAL"
onSections="SECTIONS"
def addToPositons(addOns:list,receiver:list)->list:
    for ele in np.array(addOns).flatten().tolist():
        if ele not in receiver: receiver.append( round(ele,6) )
    return list(sorted(receiver))

class Observable:
    def __init__(self):
        self._observers = []

    def add_observer(self, observer):
        """Add an observer."""
        if observer not in self._observers:
            self._observers.append(observer)

    def remove_observer(self, observer):
        """Remove an observer."""
        self._observers.remove(observer)


class Component:
    """
    Classe mere. 
    Voir CoreComponent 
    """
    def __init__(self,zMin:float,zMax:float):
        self.zMin=zMin
        self.zMax=zMax
        self._axialSubLayers_mat = [] # A quoi sa sert deja 
        self._axialAddOns_mat    = [] # A quoi sa sert deja 
        self._axialSubLayers_xs  = [] # A quoi sa sert deja 
        self._axialAddOns_xs     = [] # A quoi sa sert deja 
    def _setAxialPty_mat(self):pass
    def addToPositons(self,addOns,receiver):
        self._setAxialPty_mat()
        return addToPositons(self._axialSubLayers_mat,receiver)

class WindowHeterogeneity:
    def __init__(self,e_x,e_y,e_z,macroRegionName:str):
        self.e_x = list(sorted(e_x))
        self.e_y = list(sorted(e_y)) 
        self.e_z = list(sorted(e_z))
        self.macroRegionName = macroRegionName

    #def overlaps(self, other):return (self.e_x < other.e_x and self.e_y < other.e_y and self.e_z < other.e_z)
    def overlaps(self, other):
        """
        Check if this WindowHeterogeneity overlaps with another in 3D space.

        Parameters:
        - other (WindowHeterogeneity): The other instance to compare against.

        Returns:
        - bool: True if the two instances overlap, False otherwise.
        """
        x_overlap = not (self.e_x[1] <= other.e_x[0] or self.e_x[0] >= other.e_x[1])
        y_overlap = not (self.e_y[1] <= other.e_y[0] or self.e_y[0] >= other.e_y[1])
        z_overlap = not (self.e_z[1] <= other.e_z[0] or self.e_z[0] >= other.e_z[1])
        return x_overlap and y_overlap and z_overlap
    def is_included(self, other):
        """Check if this window is fully included in another."""
        return (
            self.e_x[0] >= other.e_x[0] and self.e_x[1] <= other.e_x[1] and
            self.e_y[0] >= other.e_y[0] and self.e_y[1] <= other.e_y[1] and
            self.e_z[0] >= other.e_z[0] and self.e_z[1] <= other.e_z[1]
        )
    def size(self):
        """Calculate the volume of the heterogeneity window."""
        return (self.e_x[1] - self.e_x[0]) * (self.e_y[1] - self.e_y[0]) * (self.e_z[1] - self.e_z[0])
    def __eq__(self, other):
        return (
            isinstance(other, WindowHeterogeneity) and
            self.e_x == other.e_x and
            self.e_y == other.e_y and
            self.e_z == other.e_z and
            self.macroRegionName == other.macroRegionName
        )

    def __hash__(self):
        return hash((tuple(self.e_x), tuple(self.e_y), tuple(self.e_z), self.macroRegionName))

    @staticmethod
    def build_hierarchy(windows):
        """
        Build a hierarchical structure of windows.

        Parameters:
        - windows (list of WindowHeterogeneity): The list of windows to process.

        Returns:
        - list of WindowHeterogeneity: A list of root windows, each with a hierarchy of children.
        """
        # Sort by size (largest first)
        windows = sorted(windows, key=lambda w: w.size(), reverse=True)

        def find_children(parent, candidates):
            """Recursively find children for a given parent window."""
            children = [w for w in candidates if w.is_included(parent) and w != parent]
            for child in children:
                remaining_candidates = [w for w in candidates if w != child]
                child.children = find_children(child, remaining_candidates)
            return children

        roots = []
        while windows:
            root = windows[0]
            windows = [w for w in windows if w != root]
            root.children = find_children(root, windows)
            roots.append(root)

        return roots

    @staticmethod
    def resolve_overlaps(windows):
        """
        Resolve overlaps among a list of WindowHeterogeneity instances by prioritizing larger windows.

        Parameters:
        - windows (list of WindowHeterogeneity): The list of windows to process.

        Returns:
        - list of WindowHeterogeneity: A filtered list of windows with resolved overlaps.
        """
        # Sort windows by size (largest first)
        windows = sorted(windows, key=lambda w: w.size(), reverse=True)

        resolved = []
        for window in windows:
            if all(not window.overlaps(w) for w in resolved):
                resolved.append(window)

        return resolve

    def __repr__(self) : return f"Xbounds: {self.e_x} , YBounds = {self.e_y} , ZBounds = {self.e_z} , xs = {self.macroRegionName}  " 


class CoreComponent(Component):
    """
    Composant d'un dessin de coeur. Utiliser pour gerer l'affectation des macro-regions si heterogeneisation necessaires. 
    Utilisable pour tout composant (ZA,Gaine,Reflecteur,tambours de controle de reactivite. 

    Philosophie: 
    -----------
     Chaque composant s'appuie sur une forme geometrique plus ou moins elementaire. L'instance de la classe generer pour cree le composant emporte avec lui toutes les
     informations necessaires pour l'insertion du composant au sein du dessin du coeur. 
     Ces informations sont generalement :
                                        - le nom du materiau 
                                        -les epaisseurs des heterogeneisations le long des axes elementaires (ex,ey,ez)
                                        -la position axiale minimale et maximale a laquelle on trouve le composant. 

    Attributs :
    ---------
    - shape   (closeGeometry) : Forme sur laquelle s'appuie le composant du coeur 
    - matName (str)           : Nom du materiau composant le materiaux 
    - ex      (float)         : Demi-epaisseur du composant selon l'axe x. 
    - ey      (float)         : Demi-epaisseur du composant selon l'axe y. 
    - ez      (float)         : Demi-epaisseur du composant selon l'axe z. 
    - zMin    (float)         : Position axiale minimale a laquelle on trouve le composant 
    - zMax    (float)         : Position axiale minimale a laquelle on trouve le composant 
    """
    def __init__(self,shape:cg,ex:float,ey:float,ez:float):
        self.shape=shape 
        self.matName=self.shape.get_data("material")
        self.ex , self.ey, self.ez = ex , ey ,ez 
        self.sectionHeterogeneity=[]
        self.zMax,self.zMin = None,None
        super().__init__(self.zMax,self.zMin)

    def setAxialPositions(self,translation):
        self.zMin=round(-self.shape.get_height()/2+translation ,6)
        self.zMax=round( self.shape.get_height()/2+translation ,6)

    def get_height(self):return self.shape.get_height()
    def get_data(self,dataName:str):return self.shape.get_data(dataName)
    def getXpos(self):return self.get_data("radius")
    def setHeterogeneity(self,e_x,e_z):
        """
        Ajoute une hetereogeneite d'epaiseur e_x selon l'axe x, et d'epaissseur selon l'axe e_z. 
        """
        self.sectionHeterogeneity.append((e_x,e_z))
    

    def getSectionsShapes(self)->cg:
        """
        Renvoie la forme geometrique ( object closeGeometry d'INCA) avec la prise en compte des heterogeneites des macro-regions. Ce resultat permettra de definir le dessins 
        contenant les affectations spatiales de macroRegions. 
        """
        heterogeneities=[]
        if self.sectionHeterogeneity==[] : 
            s= Cylinder(f"{self.shape.name}_0", self.shape.get_data("radius"),self.shape.get_height() , material = f"{self.matName}_0", xset="{self.matName}_0")
            return [s] 
        for iele,ele in enumerate(self.sectionHeterogeneity):
            ex,ez=ele
            radius,height=self.shape.get_data("radius") , self.shape.get_height()
            heterogeneities.append( Cylinder(f"{self.shape.name}_{iele}",radius-ex,height-2*ez,material=f"{self.matName}_{iele}" , xset=f"{self.matName}_{iele}") )
        s= Cylinder(f"{self.shape.name}_{iele+1}", self.shape.get_data("radius"),self.shape.get_height() , material = f"{self.matName}_{iele+1}", xset="{self.matName}_{iele+1}")
        heterogeneities.append(s)
        return heterogeneities

class MaterialLayer(Observable):
    """
    Classe permettant de gerer par couche axiale les proprietes des differents materiaux mis en jeu.  
    Attribut:
    --------
    - materialName    (str)   : Nom du materiau remplissant la couche 
    - heterogeneities (list)  : liste de couples e_x,e_y des heterogeneites sur les macroRegions
    - zMin            (float) : Position minimale axiale de la couche 
    - zMax            (float) : Position maximale axiale de la couche 
    - thickness       (float) : epaisseur de la couche 
    """
    def __init__(self, materialName:str,zMin:float,thickness:float,xmin:float,length:float,thetaMin:float,dTheta:float,symetric=False):
        Observable.__init__(self)
        #if thickness <= 0:raise ValueError("L'epaisseur doit etre une valeur positive.")
        self.materialName        = materialName
        # axial Data 
        #------------
        self.setZdata(zMin,thickness)
        # 0 axis data 
        #------------
        self.setXData(xmin,length)
        # 1 axis data 
        #-------------
        self.setYdata(thetaMin,dTheta)

        # others 
        #------
        self.symetric            = symetric 
        self.heterogeneities     = []  # List to store (e_x, e_z) pairs

        self.masterHeterogeneity = None 
        self.is_ready = False  # Track readiness internally

    def __setAxisData(self,start:float,thickness :float) :
        if start !=None and thickness != None : return start, round(start+thickness,6)
        else                                  : return start , thickness

    def setZdata(self,zMin:float     ,thickness:float)       : self.zMin     , self.zMax     = self.__setAxisData(zMin,thickness)
    def setYdata(self,thetaMin:float ,thickness:float)       : self.thetaMin , self.thetaMax = self.__setAxisData(thetaMin,thickness) 
    def setXData(self,xmin:float     ,length:float   )       : self.xmin     , self.xmax     = self.__setAxisData(xmin,length) 
    def setDataOnAxis(self, start:float , thickness : float, iaxis:int ) :
        if iaxis not in [0, 1, 2]:raise ValueError("iaxis must be 0 (x-axis), 1 (y-axis), or 2 (z-axis).")
        elif iaxis ==0  : self.setXData(start,thickness)
        elif iaxis ==1  : self.setYData(start,thickness)
        elif iaxis ==2  : self.setZData(start,thickness)
    def getAxisLength(self,iaxis):
        if iaxis==0 : return round(self.xmax-self.xmin,6)
        elif iaxis==1 : return round(self.thetaMax-self.thetaMin,6)
        elif iaxis==2 : return round(self.zMax-self.zMin,6)
        

    def _getBounds(self): return [self.xmin , self.xmax ] , [self.thetaMin,self.thetaMax] , [self.zMin,self.zMax]

    def uppdateMaterialName(self,newName:str): self.materialName = newName

    def addHeterogeneity(self,e_x:float,e_y:float,e_z:float,macroRegionIndex:int):
        """
        Add a heterogeneity to the material layer.

        Parameters:
        - e_x (float): Offset for the heterogeneity in the x-direction.
        - e_y (float): Offset for the heterogeneity in the y-direction.
        - e_z (float): Offset for the heterogeneity in the z-direction.
        - macroRegionIndex (int): Index to define the macroRegion for the heterogeneity.

        Example:
            layer = MaterialLayer("CLAD", zMin=0.0, thickness=0.5, xmin=1.0, length=2.0, thetaMin=0.0, dTheta=30.0)
            layer.addHeterogeneity(e_x=0.1, e_y=0.2, e_z=0.05, macroRegionIndex=1)
        """
        xBounds,yBounds,axialBounds = self._getBounds()
        _x,_y,_z =[ self.getAxisLength(iele)  if ele==0 else ele for iele,ele in enumerate([e_x,e_y,e_z]) ]
        #if e_x == 0 : e_x = self.getAxisLength(0)
        #if e_y==0   : e_y = self.getAxisLength(1)
        #if e_z==0   : e_z = self.getAxisLength(2) 
        heterogeneity = WindowHeterogeneity( [self.xmin,round(self.xmin+_x,6)]  
                                            ,[self.thetaMin,round(yBounds[0]+_y,6)] 
                                            ,[self.zMin,round(self.zMin+_z,6)] 
                                            ,f"{self.materialName}_{macroRegionIndex}" 
                                            )     
        if heterogeneity not in self.heterogeneities :  self.heterogeneities.append( heterogeneity ) 
        if self.symetric:
            if self.zMax - e_z <= self.zMin:raise ValueError("Symmetric heterogeneity exceeds layer bounds.")
            #if e_z !=0  : _z = -_z 
            #if e_y != 0 : _y = -_y 
            #print(_y,_z)

            heterogeneity= WindowHeterogeneity( 
                                              [self.xmin,round(self.xmin+_x,6)]
                                              ,[self.thetaMin,round(yBounds[1]-e_y,6)]
                                              ,[round(self.zMax-e_z,6),self.zMax]
                                              ,f"{self.materialName}_{macroRegionIndex}" 
                                              )     
            if heterogeneity not in self.heterogeneities :  self.heterogeneities.append( heterogeneity )  

    def getHeterogeneities(self)->list :
        return [self.masterHeterogeneity]+ self.heterogeneities
    def set_ready(self):
        """Mark the layer as ready and notify observers."""
        self.is_ready = True
        x,y,z=self._getBounds()
        self.masterHeterogeneity = WindowHeterogeneity( x,y,z, f"{self.materialName}_0") 
    def reset_ready(self):
        """Mark the layer as not ready and notify observers."""
        self.is_ready = False
        self.masterHeterogeneity=None
        self.notifyObservers(materialName=self.materialName, status="not_ready")
    def buildMPOSet(self,sectionSetter:SectionSetter):
        #for hetero in WindowHeterogeneity.resolve_overlaps(self.getHeterogeneities()):
        for hetero in WindowHeterogeneity.build_hierarchy(self.getHeterogeneities()):
            print(hetero.e_x,hetero.e_y,hetero.e_z,hetero.macroRegionName)



class Leg:
    """
    Classe permettant de gerer les branches chaudes et froides d'un systeme de boucles 

    Attributs :
    -----------

    - zMin        : position axiale minimale de la branche. Ici on s'interesse a la borne inferieur par rapport au sel.  
    - diameter    : epaisseur de la branche  
    - xmin        : position minimale selon l'axe 0 de la branche 

    Schema de principe :  
    ------------------ 
    
              |--------------------------------------------------------------           ^
              |                         GAINE                                           | e_clad       
              |-------------------------------------------------------------- <-- zMax  v 
              |                         SEL                                             ^ 
              |                                                                         | diameter 
              |                                                                         v
              |-------------------------------------------------------------- <-- zMin  ^
              |                         GAINE                                           | e_clad    
              |--------------------------------------------------------------           v
     ZONE     |
    ACTIVE    |
              |
              |

    """
    def __init__(self,zMin:float,diameter:float,xmin:float,length:float,minTheta:float,dTheta:float):
        self.zMin      = zMin 
        self.diameter  = diameter 
        self.zMax      = round(self.zMin+self.diameter,6)
        self.xmin      = xmin
        self.length    = length
        self.pipe      = MaterialLayer("SEL" ,zMin                      , self.diameter , self.xmin,self.length,minTheta,dTheta,True)
        for ele in ["cladMinus","cladPlus","upperClad","lowerClad"] : setattr(self,ele, MaterialLayer("CLAD",None,None,self.xmin,self.length,None,None)) 
        for ele in ["reflMinus","reflPlus","upperRefl","lowerRefl"] : setattr(self,ele, MaterialLayer("REFL",None,None,self.xmin,self.length,None,None)) 
        self.pipe.add_observer(self) ; self.pipe.set_ready()
        self.PIPE , self.CLAD , self.REFL = "PIPE","CLAD","REFL"
        self.UPPER , self.LOWER , self.MINUS , self.PLUS = "UPPER","LOWER","MINUS","PLUS"
        self.layerMap ={}

    def setClad(self,dTheta:float,e_clad:float):
        minTheta,maxTheta,dTheta_plugs=self.pipe.thetaMin , self.pipe.thetaMax, round(self.pipe.thetaMax-self.pipe.thetaMin+2*dTheta,6)
        for ele in ["cladPlus","cladMinus"]:getattr(self,ele).setZdata(self.zMax,round(self.diameter+e_clad,6))
        for ele in ["lowerClad","upperClad"] : getattr(self,ele).setYdata(minTheta-dTheta,dTheta_plugs)
        self.lowerClad.setZdata(self.zMin-e_clad,e_clad)
        self.upperClad.setZdata(self.zMax,e_clad)

        self.cladPlus.setYdata(maxTheta,dTheta)   # theta superieur a thetaMax du sel
        self.cladMinus.setYdata(minTheta,-dTheta) # theta inferieur a thetaMin du sel 

        for layer in [self.lowerClad,self.upperClad,self.cladMinus,self.cladPlus] : 
            layer.add_observer(self)
            layer.set_ready()
        if self.is_complete() : self.updateLayers()

    def setRefl(self , dTheta :float ,e_refl:float):
        for ele in [self.upperClad,self.lowerClad,self.cladMinus,self.cladPlus] : 
            if not ele.is_ready : raise ValueError("Veuillez renseignez les donnees sur la gaine avant de definir celles du reflecteur")

        minTheta , maxTheta ,dTheta_plugs = self.lowerClad.thetaMin, self.lowerClad.thetaMax , round(self.upperClad.thetaMax-self.upperClad.thetaMin+2*dTheta,6)

        for ele in ["reflPlus","reflMinus"]:getattr(self,ele).setZdata(self.upperClad.zMax,round(self.diameter+e_refl,6))
        for ele in ["lowerRefl","upperRefl"] : getattr(self,ele).setYdata(minTheta-dTheta,dTheta_plugs)
        self.lowerRefl.setZdata(round( self.lowerClad.zMin-e_refl,6) ,e_refl)
        self.upperRefl.setZdata( self.upperClad.zMax,e_refl)
        self.reflPlus.setYdata(   maxTheta, dTheta  )        # theta superieur a thetaMax de la gaine 
        self.reflMinus.setYdata(  minTheta,-dTheta  )        # theta inferieur a thetaMin de la gaine 
        for layer in [self.upperRefl,self.lowerRefl,self.reflMinus,self.reflPlus] : 
            layer.add_observer(self)
            layer.set_ready()
        if self.is_complete() : self.updateLayers()

        
    def updateLayers(self):
        self.layerMap ={
                         self.PIPE : self.pipe
                        ,self.CLAD :{
                                     self.LOWER : self.lowerClad 
                                    ,self.UPPER : self.upperClad 
                                    ,self.MINUS : self.cladMinus
                                    ,self.PLUS  : self.cladPlus 
                                    }
                        ,self.REFL :{ 
                                       self.LOWER :self.lowerRefl
                                      ,self.UPPER :self.upperRefl
                                      ,self.MINUS :self.reflMinus
                                      ,self.PLUS  :self.reflPlus
                                      }

                       }

    def is_complete(self):
        """Check if all clad and reflector layers are ready."""
        clad_ready = all(
            getattr(self, attr).is_ready
            for attr in ["cladPlus", "cladMinus", "upperClad", "lowerClad"]
        )
        refl_ready = all(
            getattr(self, attr).is_ready
            for attr in ["reflPlus", "reflMinus", "upperRefl", "lowerRefl"]
        )
        return clad_ready and refl_ready
    def _get_layer_by_sub(self, where: str, sub: str = None):
        """
        Retrieve the specific layer based on `where` and optionally `sub`.

        Parameters:
        - where (str): The material type (`PIPE`, `CLAD`, or `REFL`).
        - sub (str): The sub-layer type (`PLUS`, `MINUS`, `UPPER`, or `LOWER`).

        Returns:
        - MaterialLayer: The corresponding layer.
        """
        # Handle the pipe layer (no subdivisions)
        if where == self.PIPE:
            if sub is not None:raise ValueError("The 'PIPE' layer does not have subdivisions.")
            return self.pipe

        if where not in self.layerMap: raise ValueError(f"Invalid material type '{where}'. Must be one of PIPE, CLAD, or REFL.")
        if sub not in self.layerMap[where] : raise ValueError(f"Invalid sub-layer '{sub}' for material type '{where}'.")

        return self.layerMap[where][sub]

    def addHeterogeneity(self,where:str,sub:str,e_x,e_y,e_z,macroRegionIndex):
        layer=self._get_layer_by_sub(where,sub)
        layer.addHeterogeneity(e_x,e_y,e_z,macroRegionIndex)
    def buildMPOSet(self,sectionSetter):
        for layerName,layers in self.layerMap.items():
            if isinstance(layers,dict):
                for subName,sub in layers.items() : sub.buildMPOSet(sectionSetter)
            else : layers.buildMPOSet(sectionSetter)




#class SimpleLoops(CoreComponent):
class SimpleLoops:
    """
    Objet permettant de generer des boucles de manieres generiques en systemes r,theta,z.

    Philosophie :
    ------------
    L'instance de cette classe embarque avec elle :
                                                    - la longueur de la boucle jusqu'a l'echangeur.
                                                    - le nombre de boucle
                                                    - le diametre des tuyaux
                                                    - les positions axiales minimales et maximales des boucles
    A cela on considerera toujours que :
                                                    - le sel circulant dans les tuyaux est gaine
                                                    - la gaine des boucles est entoure de reflecteur.

    Attributs :
    ---------
    - nLoops    (int)   : Nombre de boucle
    - diameter  (float) : Diametre/Epaisseur des tuyaux
    - length    (float) : Longeur allant de la zone active jusqu'a la borne exterieure (selon x) de l'echangeur.

    Configuration Nominale :
    ----------------------

                            |-------------------------------------------------------------
                            |                   REFL_0
                            |
                            |-------------------------------------------------------------
                            |                   GAINE
                            |-------------------------------------------------------------
                            |                   SEL_1
                            |-------------------------------------------------------------
                            |
                            |                   SEL_0
                            |
                            |-------------------------------------------------------------
                            |                   SEL_1
                            |-------------------------------------------------------------
                            |                   GAINE
          ZONE              |-------------------------------------------------------------
        ACTIVE              |
                            |                   REFL_0
                            |-------------------------------------------------------------
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |   REFL_0      |                       REFL_1
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |               |
                            |-------------------------------------------------------------
                            |                   REFL_0
                            |
                            |-------------------------------------------------------------
                            |                   GAINE
                            |-------------------------------------------------------------
                            |                   GAINE
                            |-------------------------------------------------------------
                            |                    SEL_1
                            | -------------------------------------------------------------
                            |
                            |                    SEL_0
                            |
                            | -------------------------------------------------------------
                            |                    SEL_1
                            | -------------------------------------------------------------
                            |                    GAINE
                            |-------------------------------------------------------------
                            |
                            |                   REFL_0
                            |-------------------------------------------------------------
                            |
                            |
                            |
                            |
    Attributs :
    ----------
            - nPipes               (int)    : Nombre de tuyaux selon theta. Le nombre totale de tuyaux sera de nLoops x 2 puisque branche chaude et froide
            - diameter             (float)  : Epaisseur des tuyaux
            - length               (float)  : Longeur allant de la zone active jusqu'a la borne externe de l'echangeur
            - xmin                 (float)  : position (axe 0) minimale des boucles ( generalement rayon de la zone active)
"""
    def __init__(self,nPipes:int,length:float,diameter:float):
        #super().__init__(0,0)
        self.nPipes   = nPipes 
        self.diameter = diameter
        self.length   = length
        self.xmin     = None
        self.coldLeg  = None 
        self.hotLeg   = None 
        self.zMin,self.zMax= None, None

    def setXmin(self,xmin):self.xmin=xmin
    def setColdLeg(self,coldLeg:Leg) : 
        self.coldLeg=coldLeg
        self.zMin = self.coldLeg.zMin
    def setHotLeg(self,hotLeg:Leg):
        self.hotLeg = hotLeg
        self.zMax = self.hotLeg.zMax
    def buildMPOSet(self,sectionSetter:SectionSetter.SectionSetter,X,Y,Z):
        pass 



#class SimpleLoops(Component):
#    """
#    Objet permettant de generer des boucles de manieres generiques en systemes r,theta,z. 
#
#    Philosophie :
#    ------------
#    L'instance de cette classe embarque avec elle : 
#                                                    - la longueur de la boucle jusqu'a l'echangeur. 
#                                                    - le nombre de boucle 
#                                                    - le diametre des tuyaux 
#                                                    - les positions axiales minimales et maximales des boucles 
#    A cela on considerera toujours que :
#                                                    - le sel circulant dans les tuyaux est gaine 
#                                                    - la gaine des boucles est entoure de reflecteur. 
#
#    Attributs :
#    ---------
#    - nLoops    (int)   : Nombre de boucle 
#    - diameter  (float) : Diametre/Epaisseur des tuyaux 
#    - length    (float) : Longeur allant de la zone active jusqu'a la borne exterieure (selon x) de l'echangeur. 
#    - e_clad    (float) : Epaisseur de la macor-region associe a la gaine 
#    - e_salt    (float) : Epaisseur de la premiere macro-region de sel (a l'interface avec la gaine). 
#    - e_refl    (float) : Epaisseur de la premiere macro-region de reflecteur ( a l'interface avec la gaine)
#
#    Configuration Nominale :
#    ----------------------
#
#                            |-------------------------------------------------------------
#                            |                   REFL_0
#                            |
#                            |-------------------------------------------------------------
#                            |                   GAINE 
#                            |-------------------------------------------------------------
#                            |                   SEL_1
#                            |-------------------------------------------------------------
#                            |                   
#                            |                   SEL_0
#                            |                   
#                            |-------------------------------------------------------------
#                            |                   SEL_1
#                            |-------------------------------------------------------------
#                            |                   GAINE 
#          ZONE              |-------------------------------------------------------------
#        ACTIVE              |
#                            |                   REFL_0
#                            |-------------------------------------------------------------
#                            |               |
#                            |               |                   
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |   REFL_0      |                       REFL_1
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |               |
#                            |-------------------------------------------------------------
#                            |                   REFL_0
#                            |
#                            |-------------------------------------------------------------
#                            |                   GAINE 
#                            |-------------------------------------------------------------
#                            |                   GAINE 
#                            |-------------------------------------------------------------
#                            |                    SEL_1
#                            | -------------------------------------------------------------
#                            |                    
#                            |                    SEL_0
#                            |                    
#                            | -------------------------------------------------------------
#                            |                    SEL_1
#                            | -------------------------------------------------------------
#                            |                    GAINE 
#                            |-------------------------------------------------------------
#                            |
#                            |                   REFL_0
#                            |-------------------------------------------------------------
#                            |
#                            |
#                            |
#                            |
#    Attributs :
#    ----------
#            - nLoops               (int)   : Nombre de tuyaux selon theta. Le nombre totale de tuyaux sera de nLoops x 2 puisque branche chaude et froide 
#            - diameter             (float) : Epaisseur des tuyaux 
#            - length               (float) : Longeur allant de la zone active jusqu'a la borne externe de l'echangeur 
#            - e_clad               (float) : epaisseur de l'unique macro-region de la gaine 
#            - e_refl               (float) : epaisseur de la premiere region de reflecteur (REFL_0) a l'interface avec la gaine 
#            - e_salt               (float) : epaisseur de la premiere region de sel (SEL_1) a l'interface avec la gaine 
#            - SALT                 (str )  : Nom du materiau associe au sel dans les boucles     
#            - REFL                 (str )  : Pointeur sur le reflecteur     
#            - CLAD                 (str )  : Pointeur sur la gaine 
#            - ALL                  (str )  : T qui  
#            - saltWindow           (list)  : Contient les heterogeneites des macro-regions   
#            - cladWindow           (list)  : Contient les heterogeneites des macro-regions   
#            - reflWindow           (list)  : Contient les heterogeneites des macro-regions   
#            - xmin                (float)  : position (axe 0) minimale des boucles ( generalement rayon de la zone active) 
#    """       
#    def __init__(self,nLoops:int,length:float,diameter:float,zMin:float,zMax:float,e_salt=1.42,e_clad=1.,e_refl=5):
#        super().__init__(zMin,zMax)
#        self.nLoops=nLoops
#        self.diameter=diameter
#        self.length=length
#        self.e_clad=e_clad
#        self.e_refl=e_refl
#        self.e_salt=e_salt
#        self.SALT ="SALT"
#        self.REFL="REFL"
#        self.CLAD="CLAD"
#        self.ALL="ALL"
#        self.saltWindow=[]
#        self.cladWindow=[]
#        self.reflWindow=[]
#        self.xmin=None
#        #self.sectionHeterogeneity=[]
#    def setXmin(self,xmin):self.xmin=xmin
#
#
#    def _setAxialPty_mat(self):
#        """
#        Fournit les positions axiales d'insertion des differents elements de la boucle ( tuyaux haut et bas et gaine). 
#        Cette methode est appele pour l'affectation des materiaux ; definition du dessin du coeur. 
#        """
#        self._axialSubLayers_mat=[ [self.zMin+ele for ele in [-self.e_clad , 0 ,self.diameter , self.diameter+self.e_clad ] ]
#                                 , [self.zMax-self.diameter+ele 
#                                     for ele in [-self.e_clad , 0 ,self.diameter , self.diameter+self.e_clad ]
#                                   ]
#                                 ]
#    def noName(self,start,_slices,which):
#        """
#        Parametres :
#        -----------
#        - start   (float) : Borne de depart 
#        - _slices (list)  : liste de position axiale (qui sont soit strictement superieur ou inferieur au parametres start. Depend si branche chaude ou froide.) 
#        - which   (list)  : Attribut contenant les couples d'heterogeneites associes aux tuyaux , aux gaines ou au reflecteur. 
#        """
#        minIndex=[abs(val-start) for val in _slices].index(min([abs(val-start) for val in _slices])) # On cherche le point le plus proches de la borne de depart. 
#
#        # Le point le plus proche et la point de depart forment une borne. On doit donc renseigner la borne en respectant l'ordre. 
#        zBounds=list(sorted([round(start,6),round(_slices[minIndex],6)]))
#        #  On ajoute cette nouvelle heterogeneite a la fenetre icialement fournit en argument de la methode . 
#        which.append( [ [0,self.length] , [ round(zBounds[0],6),round(zBounds[1],6) ] ] )
#
#
#    # /!\ _setSaltHeterogeneity et _setCladHeterogeneity peuvent etre generalise . Il est peut etre judicieux de faire appel a un design pattern pour cela
#    # peut etre meme un design pattern 
#    # Un attribut branche chaude branche froide peut etre interessant ( mais que ce passe t il si on considere plusieurs entreee ? )
#    # 
#
#    def _setSaltHeterogeneity(self,e_x,e_z):
#        """
#        Ajoute une nouvelle heterogeneite dans le sel des boucles d'epaisseur e_x (selon l'axe 0) et e_z (selon l'axe 2).
#        La methode s'assure qu'il n'y as pas de doublon dans les directives d'heterogeneisation. 
#        De plus elle s'assure que les bornes soient entrees dans l'ordre croissant . 
#        /!\ Pour l'instant la methode ne fonctionne que dans l'heterogeneisation axiale des tuyaux. 
#        L'heterogenisation s'applique de maniere symetrique sur la branche chaude et froide. 
#        """
#        slices=self.getSectionHetero(self.saltWindow,2)
#        zMin_coldLeg , zMax_coldleg = self.zMin                , self.zMin+self.diameter 
#        zMin_hotLeg  , zMax_hotLeg  = self.zMax-self.diameter  , self.zMax
#                # Branche froide                                                                  |           Branche chaude 
#        #addOns=[ self.zMin,self.zMin+e_z , self.zMin-e_z + self.diameter ,self.zMin+self.diameter]+[ self.zMax-self.diameter , self.zMax-self.diameter+e_z,self.zMax-e_z,self.zMax ] 
#               #                      BRANCHE FROIDE                                    |                              BRANCHE CHAUDE  
#        addOns=[ zMin_coldLeg, zMin_coldleg + e_z , zMax_coldleg-e_z, zMax_coldleg      ,        zMin_hotLeg, zMin_hotLeg+ e_z, zMax_hotLeg-e_z,zMax_hotLeg]
#        if slices==[]: slices=addOns
#        else : slices=addToPositons(addOns,[ele for b in slices for ele in b ])
#        for ele in [self.zMax-e_z,self.zMin+self.diameter-e_z]:self.noName(ele,[val for val in slices if val>ele], self.saltWindow) # Branche chaude 
#        for ele in [self.zMin+e_z,self.zMax-self.diameter+e_z]:self.noName(ele,[val for val in slices if val<ele], self.saltWindow) # Branche froide 
#
#    def _setCladHeterogeneity(self,e_x,e_z):
#        slices=self.getSectionHetero(self.cladWindow,2)
#        zMin_coldLeg , zMax_coldleg = self.zMin-self.e_clad                , self.zMax_coldleg+self.diameter+self.e_clad
#        zMin_hotLeg  , zMax_hotLeg  = self.zMax-self.diameter-self.e_clad  , self.zMax+self.e_clad 
#        #zMin_lower,zMax_lower=self.zMin-self.e_clad              , self.zMin+self.diameter+self.e_clad
#        #zMin_upper,zMax_upper=self.zMax-self.diameter-self.e_clad, self.zMax+self.e_clad
#        #addOns=[ zMin_coldLeg ,self.zMin,self.zMax,
#        addOns=[self.Zmin,zMin_lower,zMax_lower,self.zMin+self.diameter
#               ,self.zMax,self.zMax-self.diameter,zMin_upper,zMax_upper]
#        if e_z!=0:
#            for ele in [zMin_lower-e_z,zMax_lower+e_z,zMin_upper-e_z,zMax_upper+e_z]:addOns.append(ele)
#        if slices==[]:slices=list(sorted(addOns))
#        else:slices=addToPositons(list(sorted(addOns)),[ele for b in slices for ele in b ])
#
#
#    def setHeterogeneity(self,e_x,e_z,where:str):
#        if where==self.SALT : self._setSaltHeterogeneity(e_x,e_z) 
#        if where==self.CLAD:  self._setCladHeterogeneity(e_x,e_z)
#
#
#
#    def getSectionHetero(self,whichWindow,dirNum):
#        dirPath=[0,None,1]
#        if whichWindow==[] : return whichWindow
#        else               : return [ ele[dirPath[dirNum]] for w in [whichWindow] for ele in w]
#
#    def _setSectionsOnSalt(self,zmin,zmax,sectionSetter:SectionSetter.SectionSetter):
#        _z=list(sorted( [val for ele in self.getSectionHetero(self.saltWindow,2) for val in ele if val>=zmin and val<=zmax] ))
#        mediane= round( np.median( np.array( [i for i,val in enumerate(_z)] )) )
#        inside,outside_a,outside_b=_z[mediane-1:mediane+1] , _z[:mediane] , _z[mediane:]
#        sectionSetter.setOnWindow([0,self.length],[0,360],inside,"SEL_0")
#        cnt=len(outside_a)-1
#        for outside in [outside_a,outside_b]:
#            for iele,ele in enumerate(outside[:-1]) : sectionSetter.setOnWindow([0,self.length],[0,360],[ele,outside[iele+1]],f"SEL_{cnt-iele}")
#
#
#    def buildMPOSet(self,sectionSetter:SectionSetter.SectionSetter,X,Y,Z):
#        # 1. Gestion de la boucle basse
#        self._setSectionsOnSalt(self.zMin,round(self.zMin+self.diameter,6),sectionSetter)
#        self._setSectionsOnSalt(round(self.zMax-self.diameter,6),self.zMax,sectionSetter)
#
#        




                    
                

        


