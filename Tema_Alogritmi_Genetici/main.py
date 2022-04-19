import random
import math

import numpy
import copy

f = open("input.txt","r")
g = open("output.txt","w")

#in fisierul de intrare vom avea urmatorul format
# 1 - dimensiunea populatiei
# 2 - 2 valori ce reprezinta domeniul de definitie
# 3 - 3 coeficienti ai polinomului de gradul 2
# 4 - precizia cu care se lucreaza
# 5 - probabilitatea de recombinare
# 6 - probabilitatea de mutatie
# 7 - numarul de etape al algoritmului

sir = f.readline().strip()
dimensiunePop = int(sir)
sir = f.readline().strip().split(" ")
lowerBound = float(sir[0])
upperBound = float(sir[1])
sir = f.readline().strip().split(" ")
coeficienti = [float(x) for x in sir]
precizie = int(f.readline().strip())
probRecombinare = float(f.readline().strip())
probMutatie = float(f.readline().strip())
nrEtape = int(f.readline().strip())

lungimeCrom = round(math.log((upperBound-lowerBound)*pow(10,precizie),2))

def creazaCromozom():
    cromozom = ""
    for idx in range(lungimeCrom):
        cromozom+= str(random.randint(0,1))
    return cromozom

def cromozomToX(cromozom):
    intreg = int(cromozom,2)
    x = lowerBound + intreg*pow(0.1,precizie)
    return round(x,precizie)
def functieFitness(x):
    if x<lowerBound or x>upperBound:
        return 0
    sol = x*x*coeficienti[0] + x*coeficienti[1]+coeficienti[2]
    return sol

# crearea primei populatii
g.write("Populatia initiala\n")
idx= 0
fitnessMax = 0
vectorCromozomi = [[]]
sumaF = 0
while idx != dimensiunePop:
    cr = creazaCromozom()
    x = cromozomToX(cr)
    if x >= lowerBound and x <= upperBound:
        vectorCromozomi[0].append(cr)
        g.write(str(idx+1)+": "+cr+" x = "+str(x)+" f = "+str(functieFitness(x))+"\n")
        if functieFitness(x)>fitnessMax:
            fitnessMax=functieFitness(x)
            CromozomMax = cr
        sumaF+=functieFitness(x)
        idx+=1


vectorProbSelectie = [[]]
g.write("\n\nProbabilitati selectie\n")
for idx in range(len(vectorCromozomi[0])):
    cr = vectorCromozomi[0][idx]
    x = cromozomToX(cr)
    fit = functieFitness(x)
    probSelectie = fit/sumaF
    vectorProbSelectie[0].append(probSelectie)
    g.write("cromozom   "+str(idx+1) + " probabilitate "+str(probSelectie)+"\n")

vectorProbSelectie[0].sort()
IntervalProbSelectie=[[0]]
for idx in range(len(vectorProbSelectie[0])):
    if IntervalProbSelectie[0]==[0]:
        IntervalProbSelectie[0].append(vectorProbSelectie[0][idx])
    else:
        IntervalProbSelectie[0].append(vectorProbSelectie[0][idx]+IntervalProbSelectie[0][-1])
g.write("\nIntervale probabilitati selectate\n")
for idx in range(len(IntervalProbSelectie[0])):
    if idx%5==0:
        g.write(str(IntervalProbSelectie[0][idx])+"\n")
    else:
        g.write(str(IntervalProbSelectie[0][idx]) + " ")

def cautaIndexBinar(intervaleProb,u,stanga,dreapta):
    if stanga>dreapta:
        return dreapta
    else:
        mijloc = int((stanga+dreapta)/2)
        if u == intervaleProb[mijloc]:
            return mijloc
        if u<intervaleProb[mijloc]:
            return cautaIndexBinar(intervaleProb,u,stanga,mijloc-1)
        return cautaIndexBinar(intervaleProb,u,mijloc+1,dreapta)

CromDupaSelectie = [[]]
for idx in range(dimensiunePop):
    u = random.random()
    index = cautaIndexBinar(IntervalProbSelectie[0],u,0,len(IntervalProbSelectie[0]))
    g.write("\nu = "+str(u)+" selectam cromozomul "+str(index+1))
    CromDupaSelectie[0].append(vectorCromozomi[0][index])

g.write("\nDupa selectie: \n")
for idx in range(len(CromDupaSelectie[0])):
    cr = CromDupaSelectie[0][idx]
    x = cromozomToX(cr)
    fit = functieFitness(x)
    g.write(str(idx+1)+": "+CromDupaSelectie[0][idx]+" x = "+str(x)+" f = "+str(fit)+"\n")

vectorCromIncrucisare = [[]]
g.write("\nProbabilitatea de incrucisare "+str(probRecombinare)+"\n")
for idx in range(len(CromDupaSelectie[0])):
    cr = CromDupaSelectie[0][idx]
    u = random.random()
    if u<probRecombinare:
        g.write(str(idx+1)+": "+cr+" u = "+str(u)+"<"+str(probRecombinare)+ " participa\n")
        vectorCromIncrucisare[0].append((idx,cr))
    else:
        g.write(str(idx + 1) + ": " + cr + " u = " + str(u)+"\n")

vectorCromDupaRecombinare = []
vectorCromDupaRecombinare.append(copy.deepcopy(CromDupaSelectie[0]))
while vectorCromIncrucisare[0]!=[]:
    if len(vectorCromIncrucisare[0])==1:
        vectorCromIncrucisare[0].pop()
    else:
        idx1 = random.randint(0,len(vectorCromIncrucisare[0])-1)
        idx2 = random.randint(0,len(vectorCromIncrucisare[0])-1)
        if idx2 != idx1:
            cr1=vectorCromIncrucisare[0][idx1][1]
            cr2=vectorCromIncrucisare[0][idx2][1]
            idx1_initial = vectorCromIncrucisare[0][idx1][0]
            idx2_initial = vectorCromIncrucisare[0][idx2][0]
            punctIncrucisare = random.randint(0,lungimeCrom-1)
            g.write("\nrecombinare dintre cromozomul "+str(idx1_initial+1)+" cu cromozomul "+str(idx2_initial+1)+":\n")
            g.write(cr1 +" "+cr2 +" punct "+str(punctIncrucisare)+"\nRezultat ")
            cr1 = cr1[:punctIncrucisare]+cr2[punctIncrucisare:]
            cr2 = cr2[:punctIncrucisare]+cr1[punctIncrucisare:]
            vectorCromDupaRecombinare[0][idx1_initial]=cr1
            vectorCromDupaRecombinare[0][idx2_initial]=cr2
            vectorCromIncrucisare[0].pop(idx1)
            vectorCromIncrucisare[0].pop(idx2-1)
            g.write(cr1+ "  "+ cr2+"\n")
        else:
            vectorCromIncrucisare[0].pop(idx1)

g.write("\nDupa recombinare:\n")
for idx in range(len(vectorCromDupaRecombinare[0])):
    cr = vectorCromDupaRecombinare[0][idx]
    x= cromozomToX(cr)
    fit = functieFitness(x)
    g.write(str(idx + 1) + ": " + cr + " x = " + str(x) + " f = " + str(fit) +"\n")

g.write("\nProbabilitatea de mutatie pentru fiecare gena "+str(probMutatie))
g.write("\nAu fost modificati cromozomii:\n")
cromMutatie = [[]]
for idx in range(len(vectorCromDupaRecombinare[0])):
    cr = vectorCromDupaRecombinare[0][idx]
    ok=0
    for idx2 in range(len(cr)):
        u = random.random()
        if u<probMutatie:
            if ok==0:
                cromMutatie[0].append(idx)
                ok=1
            if cr[idx2] == '0':
                cr=cr[:idx2]+'1'+cr[idx2+1:]
            else:
                cr=cr[:idx2]+'0'+cr[idx2+1:]
    vectorCromDupaRecombinare[0][idx]=cr

for idx in range(len(cromMutatie[0])):
    g.write(str(cromMutatie[0][idx]+1)+"\n")

for idx in range(len(vectorCromDupaRecombinare[0])):
    cr =vectorCromDupaRecombinare[0][idx]
    x = cromozomToX(cr)
    ok = 0
    if functieFitness(x)>fitnessMax:
        ok=1
        break
if ok ==0:
    vectorCromDupaRecombinare[0][0]=CromozomMax

def PreiaMax(populatie):
    f_max =0
    for idx in range(len(populatie)):
        cr =populatie[idx]
        x = cromozomToX(cr)
        fit = functieFitness(x)
        if fit>f_max:
            f_max=fit
    return f_max


g.write("Dupa Mutatie:\n")
for idx in range(len(vectorCromDupaRecombinare[0])):
    cr = vectorCromDupaRecombinare[0][idx]
    x = cromozomToX(cr)
    fit = functieFitness(x)
    g.write(str(idx + 1) + ": " + cr + " x = " + str(x) + " f = " + str(fit) + "\n")

g.write("\n\n\n\nEvolutia Maximului:\n")
g.write(str(PreiaMax(vectorCromDupaRecombinare[0]))+"\n")


def SumaFitnes(populatie):
    suma=0
    for idx in range(len(populatie)):
        cr =populatie[idx]
        x = cromozomToX(cr)
        fit = functieFitness(x)
        suma+=fit
    return suma

def creaza_populatie_noua(lista_crom):
    sumaF = SumaFitnes(lista_crom)
    fit_max= 0
    lista_prob_selectie = []
    for idx in range(len(lista_crom)):
        cr = lista_crom[idx]
        x = cromozomToX(cr)
        fit = functieFitness(x)
        if fit >fit_max:
            cromMax = lista_crom[idx]
            fit_max=fit
        probSelectie = fit / sumaF
        lista_prob_selectie.append(probSelectie)
    lista_prob_selectie.sort()
    interval_prob_selectie=[0]
    for idx in range(len(lista_prob_selectie)):
        if len(interval_prob_selectie) == 1:
            interval_prob_selectie.append(lista_prob_selectie[idx])
        else:
            interval_prob_selectie.append(lista_prob_selectie[idx] + interval_prob_selectie[-1])
    crom_dupa_selectie = []
    for idx in range(dimensiunePop):
        u = random.random()
        index = cautaIndexBinar(interval_prob_selectie, u, 0, len(interval_prob_selectie))
        crom_dupa_selectie.append(lista_crom[index])


    vectorCromIncrucisare = []
    for idx in range(len(crom_dupa_selectie)):
        cr = crom_dupa_selectie[idx]
        u = random.random()
        if u < probRecombinare:
            vectorCromIncrucisare.append((idx, cr))

    vectorCromDupaRecombinare=copy.deepcopy(crom_dupa_selectie)
    while vectorCromIncrucisare != []:
        if len(vectorCromIncrucisare) == 1:
            vectorCromIncrucisare.pop()
        else:
            idx1 = random.randint(0, len(vectorCromIncrucisare) - 1)
            idx2 = random.randint(0, len(vectorCromIncrucisare) - 1)
            if idx2 != idx1:
                cr1 = vectorCromIncrucisare[idx1][1]
                cr2 = vectorCromIncrucisare[idx2][1]
                idx1_initial = vectorCromIncrucisare[idx1][0]
                idx2_initial = vectorCromIncrucisare[idx2][0]
                punctIncrucisare = random.randint(0, lungimeCrom - 1)
                cr1 = cr1[:punctIncrucisare] + cr2[punctIncrucisare:]
                cr2 = cr2[:punctIncrucisare] + cr1[punctIncrucisare:]
                vectorCromDupaRecombinare[idx1_initial] = cr1
                vectorCromDupaRecombinare[idx2_initial] = cr2
                vectorCromIncrucisare.pop(idx1)
                vectorCromIncrucisare.pop(idx2 - 1)
            else:
                vectorCromIncrucisare.pop(idx1)

    for idx in range(len(vectorCromDupaRecombinare)):
        cr = vectorCromDupaRecombinare[idx]
        for idx2 in range(len(cr)):
            u = random.random()
            if u < probMutatie:
                if cr[idx2] == '0':
                    cr = cr[:idx2] + '1' + cr[idx2 + 1:]
                else:
                    cr = cr[:idx2] + '0' + cr[idx2 + 1:]
        vectorCromDupaRecombinare[idx] = cr

    if fit_max>PreiaMax(vectorCromDupaRecombinare):
        vectorCromDupaRecombinare[0] = cromMax
    return vectorCromDupaRecombinare

idx = 0
pop = vectorCromDupaRecombinare[0]
while idx!=nrEtape:

    pop = creaza_populatie_noua(pop)
    g.write(str(PreiaMax(pop))+"\n")
    idx+=1








