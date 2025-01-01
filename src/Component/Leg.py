import Conception.Berserker.src.SectionSetter as SectionSetter
import Conception.Berserker.src.Component.MaterialLayer  as MaterialLayer
import Conception.Berserker.src.Component.WindowHeterogeneity as WindowHeterogeneity 
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

