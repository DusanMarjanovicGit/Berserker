import Conception.Berserker.src.Component.Observable as Observable
import Conception.Berserker.src.ShapeSlicer as ShapeSlicer
from Conception.Berserker.src.SectionSetter import SectionSetter 
from Conception.Berserker.src.Component.WindowHeterogeneity import WindowHeterogeneity

class MaterialLayer(Observable.Observable):
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
        Observable.Observable.__init__(self)
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
        heterogeneity = WindowHeterogeneity( [self.xmin,round(self.xmin+_x,6)]
                                            ,[self.thetaMin,round(yBounds[0]+_y,6)]
                                            ,[self.zMin,round(self.zMin+_z,6)]
                                            ,f"{self.materialName}_{macroRegionIndex}"
                                            )
        if heterogeneity not in self.heterogeneities :  self.heterogeneities.append( heterogeneity )
        if self.symetric:
            if self.zMax - e_z <= self.zMin:raise ValueError("Symmetric heterogeneity exceeds layer bounds.")
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
        for hetero in WindowHeterogeneity.build_hierarchy(self.getHeterogeneities()):
            print(hetero.e_x,hetero.e_y,hetero.e_z,hetero.macroRegionName)



