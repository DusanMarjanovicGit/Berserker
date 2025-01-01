import inca 
import numpy as np
from inca.component.geometry.closegeometry        import CloseGeometry
from inca.component.geometry.tube.squaretube      import SquareTube
from inca.component.geometry.tube.cube            import Cube 
from inca.component.geometry.tube.rectangulartube import RectangularTube 
from inca.component.geometry.tube.hexagonaltube   import HexagonalTube 
from inca.component.geometry.tube.octogonaltube   import OctogonalTube 
from inca.component.geometry.tube.polygonaltube   import PolygonalTube

from inca.component.geometry.tube.cylinder        import Cylinder
from inca.component.geometry.tube.cone            import Cone 
from inca.component.geometry.tube.sphere          import Sphere 
import Physics.Space.GeometryManager              as GeometryManager
from Physics.Space.GeometryManager               import Caracteristic
import Physics.Space.slicingTool                  as slt
import Physics.Space.GeometricMesh                as GeometricMesh
import random
import transformation

class BoundingBox:
    def __init__(self,shape:CloseGeometry):
        self.computeBox(shape)

    def computeBox(self,shape):
        """
        NOn c'est pas la bonne approche. Plutot utiliser x,y,z comme methode directement dispo dans inca
        """
        _shape=shape._shape
        self.dx,self.dy,self.dz=  _shape.GetDX() , _shape.GetDY(), _shape.GetDZ() 
        bbox=RectangularTube(f"bbox-{shape.name}",self.dx*2,self.dy*2,height=self.dz*2,translation=(shape.x(),shape.y(),shape.z()),rotation=shape.get_data("rotation"))
        self.volume=bbox._shape.Capacity()
        _bbox=bbox.get_bounding_box()
        self.xmin,self.ymin,self.zmin,self.xmax,self.ymax,self.zmax=_bbox[0], _bbox[1], _bbox[2], _bbox[3], _bbox[4], _bbox[5]
        self.center=(shape.x(),shape.y(),shape.z())
    def getCenter(self):return self.center


#class Caracteristic:
#    def __init__(self,name:str,value:float,requires_translation:bool,onPositions:bool=False,onThickness:bool=False):
#        self.name=name
#        self.value=value
#        self.requires_translation = requires_translation
#        self.onPositions = onPositions
#------------------------------------------------------------------------------------------------------------------------------------------
class Shape_ (GeometricMesh.MeshOnShape):
    def __init__(self,shape):
        super().__init__()
        super().setDataFromShape(shape)
        self.shape       = shape 
        self.bbox        = BoundingBox(self.shape)
        #self._dummy = Caracteristic("dummy","ymax" ,0 , False)
        self.translation = self.shape.get_data("translation")
        self.rotation    = self.shape.get_data("rotation")
        self.arguments=[0,0,0]
        self.sliced=None

    def _checkSliceQuery(self,slices,iaxis):
        if slices[-1] > self.arguments[iaxis].value : 
            raise ValueError(f" Vous ne pouvez pas decouper plus l'axe {iaxis}-{self.arguments[iaxis].name} au dela de sa valeur limite : {self.arguments[iaxis].value}") 
        if self.arguments[iaxis].value not in slices: slices.append(self.arguments[iaxis].value)

        return slices

    def _getDivideSpec(self,iaxis,slices):
        slices=self._checkSliceQuery(slices,iaxis)
        _slices=list(reversed(slices))
        inputMap=[ [self.arguments[0].value for ele in _slices] , [self.arguments[1].value for ele in _slices] , [self.arguments[2].value for ele in _slices] ]

        for iele,ele in enumerate(_slices): inputMap[iaxis][iele] = ele 
        return self.arguments[iaxis].requires_translation ,_slices , inputMap

    def _computeTranslation(self,slices:list,iaxis:int,requires_translation):
        if requires_translation: 
            translations=[]
            translations.append(self.translation)
            for i,ival in enumerate(slices[:-1]):
                t=[0,0,0]
                t[iaxis]=-(ival/2-slices[i+1]/2)
                translations.append((t[0],t[1],t[2]))
            return list(reversed(translations))

        else : return [self.translation for ele in slices]

    def addRecursely(self,slicedShapes:list,requires_translation:bool,slices:list,iaxis:int):
        translations,rotations=self._computeTranslation(slices,iaxis,requires_translation)  , [self.rotation for ele in slicedShapes]
        if len(slicedShapes)==1 : 
            slicedShapes[0].set_data("translation",translations[0])
            self.shape.add(slicedShapes[0])
            self.sliced=self.shape
        else:
            slicedShapes=list(reversed(slicedShapes))
            for igeom , geom in enumerate(slicedShapes[:-2]):
                parent=slicedShapes[igeom+1]
                slicedShapes[igeom].set_data("translation",translations[igeom])
                parent.add(slicedShapes[igeom])
            slicedShapes[-2].set_data("translation",translations[-2]) ;   slicedShapes[-1].add(slicedShapes[-2])
            self.sliced=slicedShapes[-1]
    def divide(self,iaxis,slices:list):
        requires_translation,slices_,inputMap= self._getDivideSpec(iaxis,slices) 

        self.slicedShape , required_args ,largs =  [] , []  , []

        required_args = [inputMap[iele] for iele, ele in enumerate(self.arguments) if ele.name != "dummy"]
        
        for ele in [[required_args[iarg][iele] for iarg,arg in enumerate(required_args)] for iele,ele in enumerate(slices_)]: self.slicedShape.append( self._initPattern(ele) ) 
        translations,rotations=self._computeTranslation(slices_,iaxis,requires_translation)  , [self.rotation for ele in self.slicedShape]
        for iele,ele in enumerate(self.slicedShape):
            ele.set_data("translation",list(reversed(translations))[iele])
            ele.set_data("rotation",rotations[iele])

    def applyDivision(self,iaxis,slices:list):
        requires_translation,slices_,inputMap= self._getDivideSpec(iaxis,slices) 
        self.divide(iaxis,slices)
        self.addRecursely(self.slicedShape,requires_translation,slices_,iaxis)

    def _initPattern(*args):
        pass 



class Cylinder_(Shape_,GeometricMesh.MeshOnCylinder):
    def __init__(self,shape:Cylinder):
        Shape_.__init__(self,shape)
        self.pattern=Cylinder
        GeometricMesh.MeshOnCylinder.__init__(self)
        self.setDataFromShape(self.shape)
        self.pattern=Cylinder

    def _initPattern(*args):
        self,rmax,phi2,h=args[0],args[1][0],args[1][1],args[1][2]
        return self.pattern(str(random.randint(1,1000)),rmax,h,angles=(self.phi1.value,phi2),translation=self.translation,rotation=self.rotation,internal_radius=self.rmin.value)

class SquareTube_(Shape_,GeometricMesh.MeshOnTube):
    def __init__(self,shape):
        Shape_.__init__(self,shape)
        GeometricMesh.MeshOnTube.__init__(self)
        self.setDataFromShape(self.shape)
        self.pattern=RectangularTube
    def _initPattern(*args):
        self,width,length,height = args[0],args[1][0],args[1][1],args[1][2]
        return self.pattern(str(random.randint(1,1000)) , width, length,height,translation=self.translation , rotation= self.rotation) 

class PitchedTubes_(Shape_,GeometricMesh.MeshOnPitchedTube) : 
    def __init__(self,shape):
        Shape_.__init__(self,shape)
        GeometricMesh.MeshOnPitchedTube.__init__(self)
        if isinstance(shape,OctogonalTube) : self.pattern ,self.FrameType = OctogonalTube , GeometricMesh.OCTOGONAL
        else                               : self.pattern ,self.FrameType = HexagonalTube , GeometricMesh.HEXAGONAL 
        self.setDataFromShape(self.shape)

    def _initPattern(*args):
        self,pitch,height = args[0],args[1][0],args[1][1]
        return self.pattern(str(random.randint(1,1000)) , pitch, height , translation=self.translation, rotation= self.rotation)

    def divide(self,iaxis:int,slices:list):
        if iaxis == 1 : raise ValueError( "Pour des raisons de symetrie seuls les divisions en pitch et axiales sont permises ( respectivement axes 0 et 2)")
        super().divide(iaxis,slices)


class Sphere_(Shape_,GeometricMesh.MeshOnSphere) :
    def __init__(self,shape) :
        Shape_.__init__(self,shape)
        GeometricMesh.MeshOnSphere.__init__(self)
        self.setDataFromShape(self.shape)
        self.pattern = Sphere

    def _initPattern(*args):
        self,radius,phi2,theta2=args[0],args[1][0],args[1][1],args[1][2]
        return self.pattern(str(random.randint(1,1000)) , radius, angles=(self.phi1.value,phi2) , theta_angles=(self.theta1.value,theta2) , translation=self.translation, rotation= self.rotation) 

#class Cone_(Shape_):
#    def __init__(self,shape) :
#        super().__init__(shape) 
#        self.pattern = Cone
#        self.radius_min_inf = Caracteristic("radius_min_inf",self.shape.get_data("radius_min_inf") , False , onPositions= True ) 
#        self.radius_max_inf = Caracteristic("radius_max_inf",self.shape.get_data("radius_max_inf") , False , onPositions= True ) 
#        self.radius_min_sup = Caracteristic("radius_min_sup",self.shape.get_data("radius_min_sup") , False , onPositions= True ) 

#        self.radius_max_sup = Caracteristic("radius_max_sup",self.shape.get_data("radius_max_sup") , False , onPositions= True ) 
#        self.arguments      = [ self.radius_max_inf, 
#
#    def _initPattern(self,radius_max_inf,radius_max_sup ,height):
#        return self.pattern(str(random.randint(1,1000)) ,self.radius_min_inf.value 
#                                                        , radius_max_sup 
#                                                        , self.radius_min_inf.value
#                                                        , radius_max_sup
#                                                        , height,  translation=self.translation, rotation= self.rotation) 
#    def divide



#-------------------------------------------------------------------------------------------------------------------------------------------

class ShapeFactory:
    @staticmethod
    def create_shape(shape):
        if isinstance(shape, Cylinder):
            return Cylinder_(shape)
        elif isinstance(shape, RectangularTube) or isinstance(shape,SquareTube) or isinstance(shape,Cube):
            return SquareTube_(shape)
        elif isinstance(shape, OctogonalTube) or isinstance(shape, HexagonalTube):
            return PitchedTubes_(shape)
        elif isinstance(shape, Sphere):
            return Sphere_(shape)
        else:
            raise ValueError(f"Type de forme {type(shape)} non pris en charge.")
#class Shape_ :
#    def __init__(self,shape):
#        self.shape       = shape 
#        self.bbox        = BoundingBox(self.shape)
#        self.height      = Caracteristic("height",self.shape.get_data("height"),True,onThickness=True) 
#        self._dummy = Caracteristic("dummy", 0 , False, onPositions=True)
#        self.translation = self.shape.get_data("translation")
#        self.rotation    = self.shape.get_data("rotation")
#        self.arguments=[0,0,0]
#        self.sliced=None
#
#    def _checkSliceQuery(self,slices,iaxis):
#        if slices[-1] > self.arguments[iaxis].value : 
#            raise ValueError(f" Vous ne pouvez pas decouper plus l'axe {iaxis}-{self.arguments[iaxis].name} au dela de sa valeur limite : {self.arguments[iaxis].value}") 
#        if self.arguments[iaxis].value not in slices: slices.append(self.arguments[iaxis].value)
#
#        return slices
#
#    def _getDivideSpec(self,iaxis,slices):
#        slices=self._checkSliceQuery(slices,iaxis)
#        _slices=list(reversed(slices))
#        #if self.arguments[iaxis].onThickness : _slices=list(reversed(slices))
#        #else : _slices=list(reversed(slices))
#        inputMap=[ [self.arguments[0].value for ele in _slices] , [self.arguments[1].value for ele in _slices] , [self.arguments[2].value for ele in _slices] ]
#
#        for iele,ele in enumerate(_slices): inputMap[iaxis][iele] = ele 
#        return self.arguments[iaxis].requires_translation ,_slices , inputMap
#
#    def _computeTranslation(self,slices:list,iaxis:int,requires_translation):
#        if requires_translation: 
#            translations=[]
#            translations.append(self.translation)
#            for i,ival in enumerate(slices[:-1]):
#                t=[0,0,0]
#                t[iaxis]=-(ival/2-slices[i+1]/2)
#                translations.append((t[0],t[1],t[2]))
#            return list(reversed(translations))
#
#        else : return [self.translation for ele in slices]
#
#    def addRecursely(self,slicedShapes:list,requires_translation:bool,slices:list,iaxis:int):
#        translations,rotations=self._computeTranslation(slices,iaxis,requires_translation)  , [self.rotation for ele in slicedShapes]
#        if len(slicedShapes)==1 : 
#            slicedShapes[0].set_data("translation",translations[0])
#            self.shape.add(slicedShapes[0])
#            self.sliced=self.shape
#        else:
#            slicedShapes=list(reversed(slicedShapes))
#            for igeom , geom in enumerate(slicedShapes[:-2]):
#                parent=slicedShapes[igeom+1]
#                slicedShapes[igeom].set_data("translation",translations[igeom])
#                parent.add(slicedShapes[igeom])
#            slicedShapes[-2].set_data("translation",translations[-2]) ;   slicedShapes[-1].add(slicedShapes[-2])
#            self.sliced=slicedShapes[-1]
#    def divide(self,iaxis,slices:list):
#        requires_translation,slices_,inputMap= self._getDivideSpec(iaxis,slices) 
#
#        self.slicedShape , required_args ,largs =  [] , []  , []
#
#        required_args = [inputMap[iele] for iele, ele in enumerate(self.arguments) if ele.name != "dummy"]
#        
#        for ele in [[required_args[iarg][iele] for iarg,arg in enumerate(required_args)] for iele,ele in enumerate(slices_)]: self.slicedShape.append( self._initPattern(ele) ) 
#        translations,rotations=self._computeTranslation(slices_,iaxis,requires_translation)  , [self.rotation for ele in self.slicedShape]
#        for iele,ele in enumerate(self.slicedShape):
#            ele.set_data("translation",list(reversed(translations))[iele])
#            ele.set_data("rotation",rotations[iele])
#
#    def applyDivision(self,iaxis,slices:list):
#        requires_translation,slices_,inputMap= self._getDivideSpec(iaxis,slices) 
#        self.divide(iaxis,slices)
#        self.addRecursely(self.slicedShape,requires_translation,slices_,iaxis)
#
#
#
#    def _initPattern(*args):
#        pass 
#
#
#
#
#        
#
#
#class Cylinder_(Shape_):
#    def __init__(self,shape:Cylinder):
#        super().__init__(shape)
#        self.pattern=Cylinder
#        self.rmin = Caracteristic("rmin",self.shape.get_data("internal_radius"),False,onPositions=True)
#        self.rmax = Caracteristic("rmax",self.shape.get_data("radius"),False,onPositions=True)
#        self.phi1 = Caracteristic("phi1",self.shape.get_data("angles")[0],False,onPositions=True)
#        self.phi2 = Caracteristic("phi2",self.shape.get_data("angles")[1],False,onPositions=True)
#        self.arguments=[self.rmax,self.phi2,self.height]
#
#    def _initPattern(*args):
#        self,rmax,phi2,h=args[0],args[1][0],args[1][1],args[1][2]
#        return self.pattern(str(random.randint(1,1000)),rmax,h,angles=(self.phi1.value,phi2),translation=self.translation,rotation=self.rotation,internal_radius=self.rmin.value)
#
#class SquareTube_(Shape_):
#    def __init__(self,shape):
#        super().__init__(shape)
#        self.pattern=RectangularTube
#        self.width  = Caracteristic("width" ,self.shape.get_data("width"),True,onThickness=True)
#        self.length = Caracteristic("length",self.shape.get_data("length"),True,onThickness=True)
#        self.arguments=[self.width,self.length,self.height]
#    def _initPattern(*args):
#        self,width,length,height = args[0],args[1][0],args[1][1],args[1][2]
#        return self.pattern(str(random.randint(1,1000)) , width, length,height,translation=self.translation , rotation= self.rotation) 
#
#class PitchedTubes_(Shape_) : 
#    def __init__(self,shape):
#        super().__init__(shape)
#        if isinstance(shape,OctogonalTube) : self.pattern = OctogonalTube
#        else                               : self.pattern = HexagonalTube  
#        self.pitch  = Caracteristic("pitch", self.shape.get_data("pitch") , False, onPositions=True)
#        self.arguments=[ self.pitch,self._dummy,self.height]
#
#    def _initPattern(*args):
#        self,pitch,height = args[0],args[1][0],args[1][1]
#        return self.pattern(str(random.randint(1,1000)) , pitch, height , translation=self.translation, rotation= self.rotation)
#
#    def divide(self,iaxis:int,slices:list):
#        if iaxis == 1 : raise ValueError( "Pour des raisons de symetrie seuls les divisions en pitch et axiales sont permises ( respectivement axes 0 et 2)")
#        super().divide(iaxis,slices)
#
#
#class Sphere_(Shape_) :
#    def __init__(self,shape) :
#        super().__init__(shape)
#        self.pattern = Sphere
#        self.radius  = Caracteristic("radius" , self.shape.get_data("radius")    ,False , onPositions=True )
#        self.phi1    = Caracteristic("Phi1"   , self.shape.get_data("angles")[0] ,False , onPositions=True )
#        self.phi2    = Caracteristic("Phi2"   , self.shape.get_data("angles")[1] ,False , onPositions=True )
#        self.theta1  = Caracteristic("theta1" , self.shape.get_data("theta_angles")[0] ,False , onPositions=True )
#        self.theta2  = Caracteristic("theta2" , self.shape.get_data("theta_angles")[1] ,False , onPositions=True )
#        self.arguments=[self.radius,self.phi2,self.theta2]
#
#    def _initPattern(*args):
#        self,radius,phi2,theta2=args[0],args[1][0],args[1][1],args[1][2]
#        return self.pattern(str(random.randint(1,1000)) , radius, angles=(self.phi1.value,phi2) , theta_angles=(self.theta1.value,theta2) , translation=self.translation, rotation= self.rotation) 
#
#
