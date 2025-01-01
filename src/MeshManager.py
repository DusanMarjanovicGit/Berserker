import Conception.Berserker.src.ShapeSlicer as ShapeSlicer 
import random
import copy as cp
import Physics.Space.slicingTool as slt 
import Physics.Space.GeometricMesh as GeometricMesh
from itertools import product
import Conception.Berserker.src.transformation as transformation
from tqdm import tqdm
import pickle

class InspectionPoint:
    def __init__(value,meshPoint):
        self.value    = value
        self.meshPoint=meshPoint
        self.section=None
        self.material=None
    def setValue(self,value):self.value=value
    def setMeshPoint(self,meshPoint)  : self.meshPoint = meshPoint
    def setMaterial(self,material:str): self.material  = material
    def setSection(self,section:str)  : self.section   = section

class MeshElement:
    def __init__(self,node,name:str,meshCoordinate:tuple,meshIndex:int,meshVolume,pointIn=[]):
        self.name           = name 
        self.node           = node 
        self.meshCoordinate = meshCoordinate
        self.meshIndex      = meshIndex
        self.meshVolume     = meshVolume 
        self.material       = []
        self.sections       = []
        self.pointIn        = pointIn

    def addPointIn(self,point:tuple):
        if point not in self.pointIn                 : self.pointIn.append(point)
    def setMeshCoordinate(self,meshCoordinate:tuple) : self.meshCoordinate=meshCoordinate
    def setMeshIndex(self,meshIndex)                 : self.meshIndex=meshIndex
    def setMaterial(self,material)                   : self.material=material
    def setSections(self,sections)                   : self.sections=sections
    def setMeshVolume(self,meshVolume:float)         : self.meshVolume=meshVolume
    def setNode(self,node) : self.node=node


class Inspector:
    def __init__(self,associator,geometry=None,):
        self.associator=associator
        self.foundMaterial=None
        self.inspectionResult={}
        self.geometry=geometry
        #if geometry!=None: self.inspectOnPoints(geometry,self.points)
        
    def setAssociator(self,associator):self.associator=associator

    def inspectOnPoints(self)->dict:
        inspection=[]
        nodes=self.getNodes()
        centers,capacities=self.getCenterAndVolumes(nodes)
        #for ip,p in enumerate(self.associator.association[0]):
        #    print(self.associator.association[1][ip],self.associator.meshPoints[ip])
        for p in tqdm(self.associator.subPointsCartesian,desc='inpsection') : inspection.append(self._inspectOnPoint(nodes,centers,capacities,"",p))

        d={ele[1].get_data("material"):(ele[1],[]) for ele in inspection}
        for ele in inspection:d[ele[1].get_data("material")][1].append(ele[0])
        return d
        
    def getNodes(self):return [ j for ele in self.geometry.material_nodes().values() for j in ele]

    def getCenterAndVolumes(self,nodes):
        lCenters,capacities= [],[]
        for ele in nodes:
            s=ShapeSlicer.ShapeFactory()
            s=s.create_shape(ele)
            lCenters.append(s.bbox.center)
            #capacities.append(s.bbox.volume)
            try:capacities.append(s.shape._shape.Capacity())
            except AttributeError : capacities.append(s.bbox.volume)
        return lCenters,capacities

    def _inspectOnPoint(self,nodes,centers,capacities,name:str,point:tuple):

        lisIn,datas,found=[transformation._isIn(ele,point,centers[iele]) for iele,ele in enumerate(nodes)] ,list(zip(nodes,centers,capacities)) , []
        for ii,i in enumerate(lisIn):
            if i : found.append(datas[ii])
        try:
            lspec=[(found[i][0].get_data("radius"),found[i][0].get_data("angles")[1],found[i][0].get_data("height"),round(found[i][2],5)) for i in range(len( [ele[2] for ele in found]))]
            minCapa=[ele[2] for ele in found].index(min([ele[2] for ele in found]))
            return (point,found[minCapa][0]) 
        except ValueError : 
            return (point,None)



        



class Mesher:
    def __init__(self,shape,directives,xmin=0,ymin=0,zmin=0):
        self.directives =directives 
        self._setShape(shape,xmin,ymin,zmin)
        self.meshesElement=[]
    def _setShape(self,shape,xmin,ymin,zmin):
        ShapeFactory=ShapeSlicer.ShapeFactory()  
        self.shape = ShapeFactory.create_shape(shape)
        self.shape.setMeshSlices(self.directives,xmin=xmin,ymin=ymin,zmin=zmin)

    def getNodes(self,geom):return [ j for ele in geom.material_nodes().values() for j in ele]

    def addRecursely(self,slicedShapes:list):
        slicedShapes=list(reversed(slicedShapes))
        if len(slicedShapes)==1:
            #self.shape.shape.add(slicedShapes[0])
            return slicedShapes[0] 
        else:
            for igeom , geom in enumerate(slicedShapes[:-2]):
                parent=slicedShapes[igeom+1]
                e= slicedShapes[igeom]
                parent.add(e)
            slicedShapes[-1].add(slicedShapes[-2])
            return slicedShapes[-1]

    def _addRecursely(self,slicedShapes:list,translations:list) :
        slicedShapes=list(reversed(slicedShapes))
        translations=list(reversed(translations))
        if len(slicedShapes)==1:
            slicedShapes[0].set_data("translation",translation[0])
            self.shape.shape.add(slicedShapes[0])
            return self.shape.shape
            
        else:
            for igeom , geom in enumerate(slicedShapes[:-2]):
                parent=slicedShapes[igeom+1]
                parent.set_data("translation",translations[igeom+1])
                e= slicedShapes[igeom]
                e.set_data("translation",translations[igeom])
                parent.add(e)
            slicedShapes[-1].set_data("translation",translations[-1])
            slicedShapes[-2].set_data("translation",translations[-2])
            slicedShapes[-1].add(slicedShapes[-2])
            return slicedShapes[-1]

    def translationKept(self,slicedShape:list): return [ele.get_data("translation") for ele in slicedShape]
    def rotationKept(self,slicedShape:list)   : return [ele.get_data("rotation") for ele in slicedShape]

    def _renameNodes(self,meshedGeom):
        nodes=[ j for ele in meshedGeom.material_nodes().values() for j in ele]
        for iNode, node in enumerate(nodes):node.set_data("material",str(iNode)) 


    def applyMesh(self):
        def recursiveDivide(shape, axis_index):
            if axis_index >= len(self.directives):
                # cas standard : plus de decoupage 
                return shape
            # Divise la forme selon l'axe courant 
            ShapeFactory = ShapeSlicer.ShapeFactory()
            _shape=ShapeFactory.create_shape(shape)
            indexesHandler=[]
            for iele,ele in enumerate(_shape.arguments):
                if ele.name!="dummy":indexesHandler.append(iele)
            _shape.divide(indexesHandler[axis_index], self.directives[axis_index][1:])
    
            # Translations conserve en memoire 
            translations = self.translationKept(_shape.slicedShape)
            sliced_shapes = []
    
            # Traite chaque subdivision 
            for divided_part in _shape.slicedShape:
                # Decoupage recursif selon les axes suivants  
                processed_child = recursiveDivide(divided_part, axis_index + 1)
                sliced_shapes.append(processed_child)
            # Translation remise a jour 
            for idx, part in enumerate(sliced_shapes):
                part.set_data("translation", translations[idx])
            # Ajout recursif 
            return self.addRecursely([ele for ele in sliced_shapes])
    
        # Initialisation 
        final_shape = recursiveDivide(self.shape.shape, 0)
    
        return final_shape

    def _createInspectionPoints(self,nx,ny,nz):
        association=self.shape.ROOT_Conformization(randomRange=[nx,ny,nz])
        association.verifySampling()
        return association

    def _getVolumeOnMeshes(self,inspectionPoints:GeometricMesh.MeshContainer,meshedGeom):
        inspector=Inspector(inspectionPoints,meshedGeom)
        d=inspector.inspectOnPoints()
        for mesh in self.meshesElement:
            n=len(d[mesh.node.get_data("material")][1])/len(inspectionPoints.subPointsCartesian)
            mesh.setVolume(n*self.shape.shape._shape.Capacity())

    def associateMeshAndNodes(self,inspectionPoints:GeometricMesh.MeshContainer,meshedGeom,computeVolume:bool=True):
        association = self.shape.ROOT_Conformization()
        self._renameNodes(meshedGeom)
        inspector=Inspector(association,meshedGeom)
        d=inspector.inspectOnPoints()
        for mat,val in tqdm(d.items(),desc="association-noeuds et points") : 
            meshElement=MeshElement(val[0],None,None,None,None,val[1])
            meshElement.setMeshIndex(inspector.associator.initPoints.index( inspector.associator.getInitialPoint(val[1][0])) )             # Index de ce point
            meshElement.setMeshCoordinate(inspector.associator.meshPoints[meshElement.meshIndex]                             )             # Point du maillage en coordonne dans le systeme de la forme
            self.meshesElement.append(meshElement)
        if computeVolume:self._getVolumeOnMeshes(inspectionPoints,meshedGeom)



        
        

class Linker:
    def __init__(self):
        self.geometry=None
        self.mpoGeometry=None
        self.meshedGeometry=None
        self.mesher=None 
        self.inspectionSlices=[]
        self.inspectionPoint=[]
        self.containerOnPoints=None # MeshContainer on inspectionPoints
    def setGeometry(self,geom):self.geometry=geom
    def setMpoGeom(self,geom):self.mpoGeometry=geom
    def setInspectionPoints(self,nx,ny,nz):
        self.inspectionSlices=[nx,ny,nz]

    def _createInspectionPoints(self):
        self.containerOnPoints =self.mesher._createInspectionPoints(self.inspectionSlices[0],self.inspectionSlices[1],self.inspectionSlices[2])
        for point in points:self.inspectionPoint.append( InspectionPoint(point) )

    def setMesh(self,directives,xmin=0,ymin=0,zmin=0):
        if self.geometry==None or self.mpoGeometry==None or self.inspectionSlices==[]: raise ValueError
        else:
            self.mesher=Mesher(self.geometry,directives,xmin,ymin,zmin)
            self.meshedGeometry=self.mesher.applyMesh()
            self._createInspectionPoints()
            self.mesher.associateMeshAndNodes(self.containerOnPoints,self.meshedGeometry)
    def _inspect(self):
        pass 
    def link(self):pass
    



        

