import ShapeSlicer
import random
import copy as cp 
import Physics.Space.slicingTool as slt 
import Physics.Space.GeometricMesh as GeometricMesh 
import transformation 
from tqdm import tqdm
import MeshManager
import pickle
class materialSetter:
    def __init__(self,topShape,lDirectives:list):
        self.topShape=topShape
        self.directives= lDirectives 
        self.meshedGeom=None
    def applyMesh(self):
        self.mesher=MeshManager.Mesher(self.topShape,self.directives,0,0,0)
        self.meshedGeom=self.mesher.applyMesh()
        self.mesher.associateMeshAndNodes([],self.meshedGeom,computeVolume=False)
    def _changeMeshNodePpty(self,meshNode,newValue):
        if not isinstance(newValue,str):
            newValue=newValue.get_data("material")
        meshNode.set_data("material",newValue)

    def setOnShape(self,shape):
        for meshElement in self.mesher.meshesElement:
            meshNode=meshElement.node
            s=ShapeSlicer.ShapeFactory()
            s=s.create_shape(meshNode)
            for ele in meshElement.pointIn:
                #if transformation._isIn(shape,ele,(0,0,0)):meshNode.set_data("material",shape.get_data("material"))
                if transformation._isIn(shape,ele,(0,0,0)):self._changeMeshNodePpty(meshNode,shape.get_data("material"))#meshNode.set_data("material",shape.get_data("material"))

    def setOnWindow(self,X:list,Y:list,Z:list,xsName:str):
        x_min,x_max,y_min,y_max,z_min,z_max=X[0],X[-1],Y[0],Y[-1],Z[0],Z[-1]
        for meshElement in self.mesher.meshesElement:
            _x,_y,_z=meshElement.meshCoordinate
            #if _z == z_max : print(_x,_y,_z)
            xcond=(_x>=x_min and _x<x_max)
            ycond=(_y>=y_min and _y<y_max)
            zcond=(_z>=z_min and _z<z_max)
            if xcond and ycond and zcond : 
                #print(_x,_y,_z)
                #meshElement.node.set_data("material",xsName)
                self._changeMeshNodePpty(meshElement.node,xsName)


class SectionSetter(materialSetter):
    def __init__(self,topShape,lDirectives:list):
        super().__init__(topShape,lDirectives)
    def _changeMeshNodePpty(self,meshNode,newValue):
        if not isinstance(newValue,str):
            newValue=newValue.get_data("material")
        meshNode.set_data("material",newValue)


