import inca 
import inca.component.geometry.closegeometry as cg 
from inca.component.geometry.tube.cylinder import Cylinder 
import sys 
sys.path.append("/home/dm266236/AP3/StaticCase/Coeur/Minaret/3D/Robustesse/General/ConceptionPackage")
import MeshManager
import SectionSetter
import numpy  as np
import ShapeSlicer
import Component

onMaterial="MATERIAL"
onSections="SECTIONS"
def addToPositons(addOns:list,receiver:list)->list:
    for ele in np.array(addOns).flatten().tolist():
        if ele not in receiver: receiver.append( round(ele,6) )
    return list(sorted(receiver))

class SimpleCoreCreator:
    def __init__(self,rSalt:float,H_D:float,e_clad=1.,e_refl=90,e_prot=10):

        self.fuelName,self.cladName,self.reflName,self.protName="SEL","GAINE","REFL","PROT"
        self.fuel=CoreComponent(Cylinder("0",rSalt                               , H_D*rSalt*2                                       ,material=self.fuelName ),0     ,0, 0     )
        self.clad=CoreComponent(Cylinder("1",rSalt+e_clad                        , round(self.fuel.get_height() +  (2*e_clad) , 6)   ,material=self.cladName ),e_clad,0, e_clad)
        self.refl=CoreComponent(Cylinder("2",self.clad.get_data("radius")+e_refl , round(self.clad.get_height() +  (2*e_refl) , 6)   ,material=self.reflName ),e_refl,0, e_refl)
        self.prot=CoreComponent(Cylinder("3",self.refl.get_data("radius")+e_prot , round(self.refl.get_height() +  (2*e_prot) , 6)   ,material=self.protName ),e_prot,0, e_prot)
        self._setPositions()
        self.geom=None
        self.mpoSet=None

    def resetMatName(self,oldName,newName:str):
        attNames=["fuelName","cladName","reflName","protName"]
        setattr(self,attNames[ [self.fuelName,self.cladName,self.reflName,self.protName].index(oldName) ] ,newName)

    def _setPositions(self):
        _Z=[]
        for ele in [self.fuel,self.clad,self.refl,self.prot]:_Z.append(-ele.get_height()/2)
        zmin=abs(min( _Z))
        for ele in [self.fuel,self.clad,self.refl,self.prot] : ele.setAxialPositions(zmin)

    def buildGeom(self):
        self.clad.shape.add(self.fuel.shape) ; self.refl.shape.add(self.clad.shape) ; self.prot.shape.add(self.refl.shape)
        self.geom=self.prot.shape
    
    def buildMPOSet(self):
        self.fuel.setHeterogeneity(1.21,1.42)
        self.refl.setHeterogeneity(85,85)
        nodes= [ s   for ele in [self.fuel,self.clad,self.refl,self.prot] for s in ele.getSectionsShapes()]      
        for iele,ele in enumerate(nodes[:-1]):nodes[iele+1].add(ele)
        self.mpoSet=nodes[-1]

class CoreCreator(SimpleCoreCreator):

    def __init__(self,rSalt:float,H_D:float,e_clad=1.,e_refl=90.,e_prot=10):
        super().__init__(rSalt,H_D,e_clad,e_refl,e_prot)
        self.loops=None
        self.rods=[]

    def addLoops(self,loops:SimpleLoops):
        self.loops=loops
        self.loops.xmin=self.fuel.getXpos()
    def _getMeshAssignment(self):
        X=list(sorted([round(ele.getXpos(),6) for ele in [self.fuel,self.clad,self.refl,self.prot] ]))
        Y=[45*i for i in range(1,9)]
        #Z=[ round(ele,6) for ele in [self.zMinRefl,self.zMinGaine,self.zMinSalt,self.zMaxSalt,self.zMaxGaine,self.zMaxRefl,self.zMaxProt] ]
        self._setPositions()
        Z= list(sorted( [getattr(ele,attName) for ele in [self.fuel,self.clad,self.refl,self.prot] for attName in ["zMin","zMax"] ] ))
        if self.loops !=None : Z=self.loops.addToPositons(self.loops._axialSubLayers_mat,Z)
        return X,Y,Z



    def buildGeom(self):
        if self.loops==None and self.rods==[] : super().buildGeom()
        else:
            X,Y,Z=self._getMeshAssignment()
             
            matSet=SectionSetter.materialSetter( Cylinder("PaulSayne",self.loops.length    , height=self.prot.get_height() ,material="VIDE"    ),[X,Y,Z])
            matSet.applyMesh()
            for ele in [self.prot,self.refl,self.clad,self.fuel] : matSet.setOnShape(ele.shape)
            self.loops.addLoops(matSet,onMaterial,{"SEL":[0,self.refl.getXpos()] , "ALL": [self.fuel.getXpos(),self.refl.getXpos()]})
            matSet.setOnWindow([self.fuel.getXpos() , self.clad.getXpos()] ,[0,360] , [self.clad.zMin                       , round(self.clad.zMin+self.clad.ez,6) ] ,self.clad.matName )
            matSet.setOnWindow([self.fuel.getXpos() , self.clad.getXpos()] ,[0,360] , [round(self.clad.zMax-self.clad.ez,6) , self.clad.zMax                       ] ,self.clad.matName )
            #try: self.loops.addLoops(matSet,onMaterial,{"SEL":[0,self.refl.getXpos()] , "ALL": [self.clad.getXpos(),self.refl.getXpos()]})
            #except AttributeError:pass 
            self.geom=matSet.meshedGeom

    def buildMPOSet(self):
        if (self.loops, self.rods)==(None,None):super().buildMPOSet()
        else :
            self.fuel.setHeterogeneity(1.21,1.42)
            self.refl.setHeterogeneity(85,85)
            X,Y,Z = self._getMeshAssignment()
            X=addToPositons([ele.getXpos()-ex[0] for ele in [self.prot,self.refl,self.clad,self.fuel] for ex in ele.sectionHeterogeneity],X)
            Z=addToPositons([val for ele in [self.prot,self.refl,self.clad,self.fuel] for ez in ele.sectionHeterogeneity for val in [getattr(ele,"zMin")+ez[1],getattr(ele,"zMax")-ez[1]] ],Z)
            Z=[round(ele,6) for ele in Z]
            Z=addToPositons(self.loops.getSectionHetero(self.loops.saltWindow,2),Z)
            sectionSetter=SectionSetter.SectionSetter( Cylinder("Robespierro",self.loops.length    , height=self.prot.get_height() ,material="VIDE"    ),[X,Y,Z])
            sectionSetter.applyMesh()
            nodes=[ s   for ele in [self.prot,self.refl,self.clad,self.fuel] for s in ele.getSectionsShapes()]      
            for node in nodes :sectionSetter.setOnShape(node)
            self.loops.buildMPOSet(sectionSetter,X,Y,Z)
            self.mpoSet=sectionSetter.meshedGeom
            



#    def buildMPOSet(self):
#        if (self.loops, self.rods)==(None,None):super().buildMPOSet()
#        else :
#            X,Y,Z = self._getMeshAssignment()
#            eSalt_x_mpo,egaine_mpo,eRefl_mpo,eProt_mpo = -1.21 ,0  , 5 , 0
#            eSalt_z_mpo=1.42
#            X=addToPositons([ self.rSalt-eSalt_x_mpo , self.rSalt+egaine_mpo , self.rGaine+eRefl_mpo ,self.rRefl+eProt_mpo],X)
#            Z=addToPositons([ self.zMinSalt-eSalt_z_mpo , self.zMaxGaine-self.e_clad-eSalt_z_mpo], Z)
#
#            sectionSetter = SectionSetter.SectionSetter(Cylinder("0",self.loops.length    , height=self.hProt         ,material="VIDE" , xset="VIDE"   ),[X,Y,Z])
#            sectionSetter.applyMesh()
#            #mat=["PROT","REFL","GAINE","SEL"]# Utiliser le module d'inspection pour nommer directement en fait .... 
#            specs={
#                    "PROT":[[(self.rRefl,self.zMinProt)        , (self.e_prot,self.e_prot)      ]  # 
#                            ,[(self.rRefl,self.zMaxRefl)       , (self.e_prot,self.e_prot)      ]] # 
#                           
#
#                    ,"REFL":[[(self.rGaine,self.zMinGaine)                , (eRefl_mpo,-eRefl_mpo)     ]  # REMPLACER PAR UN OBJET DE TYPE COMPOSANT FINALEMENT ...  
#                            ,[(self.rGaine,self.zMaxGaine)                , (eRefl_mpo,eRefl_mpo)      ]] #  
#                           
#
#                    ,"GAINE":[[(self.rSalt,self.zMinGaine) , (egaine_mpo,egaine_mpo)         ]  # 
#                             ,[(self.rSalt,self.zMaxSalt)  , (egaine_mpo,egaine_mpo)           ]] # 
#                            
#
#                    ,"SEL":[[(self.rSalt,self.zMinSalt)    , (eSalt_x_mpo,eSalt_z_mpo)        ]
#                            ,[(self.rSalt,self.zMaxSalt)   , (eSalt_x_mpo,eSalt_z_mpo)        ]] # 
#                  
#                           
#                 }
#            for mat,bounds in specs.items():
#                cnt=0
#                for layer in bounds:
#                    x,z=layer[0]
#                    #if (0,0) not in layer[1:] : layer=[layer[0]]+[(0,0)]+layer[1:]
#                    for sub in layer[1:]:
#                        dx,dz=sub[0],sub[1]
#                        print(list(sorted([x,x+dx])),list(sorted([z,z+dz])),f"{mat}_{cnt}")
#                        #if (dx,dz)!=(0,0):cnt+=1
#
#
#                    #for idata,data in enumerate(layer):
#                    #    x,z=data[0]
#                    #    if (0,0) not in data[1:] : data=[data[0]]+[(0,0)]+data[1:]
#                    #    for sub in data[1]:
#                    #        dx,dz=sub[0],sub[1]
#                    #        print(x+dx,z+dz,f"{mat}_{cnt}")
#                    #        if (dx,dz)!=(0,0):cnt+=1
#
#
#            #for iele,ele in enumerate(specs):
#            #    x,z = ele[0] 
#            #    cnt=0
#            #    print(ele)
#            #    if (0,0) not in ele[1:] :ele=[ele[0]]+[(0,0)]+ele[1:]
#            #    print(ele)
#            #    for sub in ele[1:]:
#            #        print(sub)
#            #        dx,dz=sub[0],sub[1]
#            #        print(mat[iele],cnt,x+dx,z+dz)
#            #        sectionSetter.setOnShape( Cylinder(str(iele+1), x+dx,height=z+dz,material=f"{mat[iele]}_{cnt}"  ) )
#            #        if (dx,dz)!=(0,0) : cnt+=1
#




                    
                

        


#class Leg(Component) :
#    """
#    Permet de definir une branche des boucles (chaudes ou froides). 
#
#    Philiosophie : 
#    -------------
#    """
#    def __init__(self, zMin,diameter,e_clad):
#        self.zMin     = zMin
#        self.zMax     = None 
#        self.diameter = diameter 
#        self.e_clad   = e_clad 
#    def _setZMax(self):
#        self.zMax = self.zMin+self.diameter+self.e_clad


         

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

