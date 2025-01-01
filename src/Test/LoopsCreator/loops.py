#import Conception.Berserker.src.Component.Leg as Leg
from Conception.Berserker.src.Component.Leg import Leg 


hotLeg  = Leg(284.84,10,46.71,110,0,20)
coldLeg = Leg(100,10,46.71,110,0,20) 
for leg in [hotLeg,coldLeg]:
    leg.setClad(2,1)
    leg.setRefl(5,5)
hotLeg.addHeterogeneity(hotLeg.PIPE,None,0,0,1.42,1)
hotLeg.addHeterogeneity(hotLeg.PIPE,None,0,2 ,0,1)
hotLeg.addHeterogeneity(hotLeg.REFL,hotLeg.UPPER,0,0,5,1)
hotLeg.addHeterogeneity(hotLeg.REFL,hotLeg.LOWER,0,0,5,1)
hotLeg.addHeterogeneity(hotLeg.REFL,hotLeg.MINUS,0,-5,0,1)
hotLeg.addHeterogeneity(hotLeg.REFL,hotLeg.PLUS,0,5,0,1)

coldLeg.addHeterogeneity(coldLeg.PIPE,None,0,0,1.42,1)
coldLeg.addHeterogeneity(coldLeg.PIPE,None,0,2 ,0,1)
coldLeg.addHeterogeneity(coldLeg.REFL,coldLeg.UPPER,0,0,5,1)
coldLeg.addHeterogeneity(coldLeg.REFL,coldLeg.LOWER,0,0,5,1)
coldLeg.addHeterogeneity(coldLeg.REFL,coldLeg.MINUS,0,-5,0,1)
coldLeg.addHeterogeneity(coldLeg.REFL,coldLeg.PLUS,0,5,0,1)


#print(hotLeg.pipe.getHeterogeneities())
#print(hotLeg.upperClad.buildMPOSet(None))
hotLeg.buildMPOSet(None)

