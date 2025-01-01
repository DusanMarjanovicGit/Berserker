from math import pi,cos,sin
import inca.component.geometry.geometry as incaGeom
import math
import ROOT
ROOT.gInterpreter.ProcessLine('#include "/home/dm266236/ToolsBox/Geometry/Homogenizer/cpp/contain.h"')
from math import cos, sin, acos, degrees, radians, atan2, pi
import numpy as np



def getEulerAngles(X, Y, Z):
    """
    Translate the rotation according (X,Y,Z) axes into Euler angles (alpha,beta,gamma).
    The choice is made to apply successively (X,Y,Z) angles in that order.
    The Euler angles definition used by ROOT is :
    alpha is the rotation angle about Z axis and is done first,
    beta is the rotation about new Y and is done second,
    gamma is the rotation angle about new Z and is done third.
    All angles are in degrees.

        :param X: First rotation according the X axe in degrees.
        :type X: double
        :param Y: Second rotation according the Y axe in degrees.
        :param Y: double
        :type Z: Third rotation according the Z axe in degrees.
        :type Z: double
        :return: Euler angles in degrees
        :rtype: tuple(double,double,double)
    """
    if Z == 0 and Y == 0 and X == 0:
        return X, Y, Z
    # Convert the angles into radians
    X = radians(X)
    Y = radians(Y)
    Z = radians(Z)
    # Define the rotation matrices
    Rx = np.array(
        [[1.0, 0.0, 0.0], [0.0, cos(X), -sin(X)], [0.0, sin(X), cos(X)]]
    )
    Ry = np.array(
        [[cos(Y), 0.0, sin(Y)], [0.0, 1.0, 0.0], [-sin(Y), 0.0, cos(Y)]]
    )
    Rz = np.array(
        [[cos(Z), -sin(Z), 0.0], [sin(Z), cos(Z), 0.0], [0.0, 0.0, 1.0]]
    )
    # Define the axes array
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])
    # Rotate the x,y,z axes
    z1, z2, z3 = Rz.dot(Ry.dot(Rx.dot(z)))
    y1, y2, y3 = Rz.dot(Ry.dot(Rx.dot(y)))
    x1, x2, x3 = Rz.dot(Ry.dot(Rx.dot(x)))
    # Calculate the Euler angles
    alpha = degrees(atan2(y3, -x3)) if z3 ** 2 != 1 else 0.0
    beta = degrees(acos(z3))
    gamma = (
        degrees(atan2(z2, z1))
        if z3 ** 2 != 1
        else degrees(z3 * atan2(x2, y2))
    )

    if alpha is None or gamma is None:
        raise ValueError

    return alpha, beta, gamma



def setRotation(rx:float,ry:float,rz:float)->list:
    """
    Retourne la matrice de rotation, en utilisant des angles d'Euler, pour connaître les coordonnés d'un point dans la nouvelles base. La matrice est retournée sous forme de liste. 

    Paamètres :
    ----------
    -rx (float) : rotation autour de l'axe x 
    -ry (float) : rotation autour de l'axe y 
    -rz (float) : rotation autour de l'axe z 

    Return :
    ------
    -rotationMatrix (list) : Matrice de rotation 
    """

    #Double_t degrad = TMath::Pi() / 180.;
    #Double_t sinphi = TMath::Sin(degrad * phi);
    #Double_t cosphi = TMath::Cos(degrad * phi);
    #Double_t sinthe = TMath::Sin(degrad * theta);
    #Double_t costhe = TMath::Cos(degrad * theta);
    #Double_t sinpsi = TMath::Sin(degrad * psi);
    #Double_t cospsi = TMath::Cos(degrad * psi);
    #fRotationMatrix[0] = cospsi * cosphi - costhe * sinphi * sinpsi;
    #fRotationMatrix[1] = -sinpsi * cosphi - costhe * sinphi * cospsi;
    #fRotationMatrix[2] = sinthe * sinphi;
    #fRotationMatrix[3] = cospsi * sinphi + costhe * cosphi * sinpsi;
    #fRotationMatrix[4] = -sinpsi * sinphi + costhe * cosphi * cospsi;
    #fRotationMatrix[5] = -sinthe * cosphi;
    #fRotationMatrix[6] = sinpsi * sinthe;
    #fRotationMatrix[7] = cospsi * sinthe;
    #fRotationMatrix[8] = costhe;

    alpha,beta,gamma=getEulerAngles(rx,ry,rz)
    #alpha, beta , gamma = rx,ry,rz
    degrad    = pi/180
    sinAlpha  = sin(degrad*alpha)
    cosAlpha  = cos(degrad*alpha)
    sinBeta   = sin(degrad*beta)
    cosBeta   = cos(degrad*beta)
    sinGamma  = sin(degrad*gamma)
    cosGamma  = cos(degrad*gamma)
    rotationMatrix    = [0]*9
    # Matrice directe 
    #rotationMatrix[0] = cosGamma * cosAlpha - cosBeta * sinAlpha * sinGamma
    #rotationMatrix[1] = -sinGamma * cosAlpha - cosBeta * sinAlpha * cosGamma
    #rotationMatrix[2] = sinBeta * sinAlpha
    #rotationMatrix[3] = cosGamma * sinAlpha + cosBeta * cosAlpha * sinGamma
    #rotationMatrix[4] = -sinGamma * sinAlpha + cosBeta * cosAlpha * cosGamma
    #rotationMatrix[5] = -sinBeta * cosAlpha
    #rotationMatrix[6] = sinGamma * sinBeta
    #rotationMatrix[7] = cosGamma * sinBeta
    #rotationMatrix[8] = cosBeta

    # Matrice directe 
    rotationMatrix[0] = cosGamma * cosAlpha - cosBeta * sinAlpha * sinGamma
    rotationMatrix[3] = -sinGamma * cosAlpha - cosBeta * sinAlpha * cosGamma
    rotationMatrix[6] = sinBeta * sinAlpha
    rotationMatrix[1] = cosGamma * sinAlpha + cosBeta * cosAlpha * sinGamma
    rotationMatrix[4] = -sinGamma * sinAlpha + cosBeta * cosAlpha * cosGamma
    rotationMatrix[7] = -sinBeta * cosAlpha
    rotationMatrix[2] = sinGamma * sinBeta
    rotationMatrix[5] = cosGamma * sinBeta
    rotationMatrix[8] = cosBeta
    return rotationMatrix
    #return rotationMatrix

def applyTranslation(translationMatrix:tuple,point:tuple)->tuple :
    return point[0]+translationMatrix[0],point[1]+translationMatrix[1],point[2]+translationMatrix[2]

def applyRotation(rotationAngles:tuple,point)->tuple:
    rotationMatrix=setRotation(rotationAngles[0],rotationAngles[1],rotationAngles[2])
    x,y,z=point
    newX= rotationMatrix[0]*x + rotationMatrix[3]*y + rotationMatrix[6]*z 
    newY= rotationMatrix[1]*x + rotationMatrix[4]*y + rotationMatrix[7]*z 
    newZ= rotationMatrix[2]*x + rotationMatrix[5]*y + rotationMatrix[8]*z 

    return newX,newY,newZ

def applyTransformation(translationMatrix:tuple, rotationAngles:tuple,point:tuple)->tuple:
    p=applyRotation(rotationAngles,point)
    return applyTranslation(translationMatrix,p)
    #p=applyTranslation(translationMatrix,point)
    #return applyRotation(rotationAngles,p)
    




def isIn(node,point:tuple) :
    # 1-- translate point
    # ------------------
    ox,oy,oz=node.get_data("translation")
    x,y,z=( point[0]-ox ,point[1]-oy ,point[2]-oz )
    #x,y,z=( point[0] ,point[1] ,point[2] )
    rx,ry,rz=node.get_data("rotation")

    rotationMatrix = setRotation(rx,ry,rz) 
    newX= rotationMatrix[0]*x + rotationMatrix[1]*y + rotationMatrix[2]*z 
    newY= rotationMatrix[3]*x + rotationMatrix[4]*y + rotationMatrix[5]*z 
    newZ= rotationMatrix[6]*x + rotationMatrix[7]*y + rotationMatrix[8]*z 
    x,y,z = (newX , newY, newZ )  
    isIn = ROOT.isIn(node._shape,x,y,z) 
    if isIn ==1 : return True 
    else : return False 
    
def _isIn(node,point:tuple,translation) :
    # 1-- translate point
    # ------------------
    ox,oy,oz=translation
    x,y,z=( point[0]-ox ,point[1]-oy ,point[2]-oz )
    #x,y,z=( point[0] ,point[1] ,point[2] )
    rx,ry,rz=node.get_data("rotation")

    rotationMatrix = setRotation(rx,ry,rz) 
    newX= rotationMatrix[0]*x + rotationMatrix[1]*y + rotationMatrix[2]*z 
    newY= rotationMatrix[3]*x + rotationMatrix[4]*y + rotationMatrix[5]*z 
    newZ= rotationMatrix[6]*x + rotationMatrix[7]*y + rotationMatrix[8]*z 
    x,y,z = (newX , newY, newZ )  
    isIn = ROOT.isIn(node._shape,x,y,z) 
    if isIn ==1 : return True 
    else : return False 
    
    










#def isIn(node,point:tuple) :
#    # 1-- translate point
#    # ------------------
#    ox,oy,oz=node.get_data("translation")
#    #print(ox,oy,oz)
#    x,y,z=point
#    rx,ry,rz=node.get_data("rotation")
#    #print(rx,ry,rz)
#
#    rotationMatrix = setRotation(rx,ry,rz) 
#    newX= rotationMatrix[0]*x + rotationMatrix[1]*y + rotationMatrix[2]*z 
#    newY= rotationMatrix[3]*x + rotationMatrix[4]*y + rotationMatrix[5]*z 
#    newZ= rotationMatrix[6]*x + rotationMatrix[7]*y + rotationMatrix[8]*z 
#    x,y,z = (newX-ox , newY-oy, newZ-oz )  
#    isIn = ROOT.isIn(node._shape,x,y,z) 
#    if isIn ==1 : return True 
#    else : return False 
