import tkinter as tk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  
from matplotlib.figure import Figure  
import re
import os
import gzip
import shutil

                                    
def liste_cds(nom_fich):
        li=[]
        f=open(nom_fich, "r")
        r= re.compile("^CDS(\t)*with_protein")
        
        for line in f:
            line=line.rstrip('\n')
            res= r.search(line)
            
            if res:
                line=line.split()
               
                li.append(line)
              
        return li       
    
cds1=liste_cds('GCA_000009985.1_ASM998v1_feature_table.txt') 
cds2=liste_cds('GCA_000014865.1_ASM1486v1_feature_table.txt')





class fenetres():
    def __init__(self, blast_file, test, cds1, cds2):
        self.blast_file=blast_file
        self.cds1=cds1
        self.cds2=cds2
        self.f = tk.Tk()
        self.f.geometry("350x200")
        self.test=test
        if self.test=='e_value':
            self.f.title("e_value ")
            champ_label = tk.Label(self.f, text="Entrez e_value maximal")
            champ_label.pack()
        if self.test=="pourcentage d'id":
            self.f.title("Pourcentage d'identité ")
            champ_label = tk.Label(self.f, text="Entrez le pourcentage d'identité minimal")
            champ_label.pack()
        if self.test=='couverture':
            self.f.title(" Taux de Couverture ")
            champ_label = tk.Label(self.f, text="Entrez taux de couverture minimal")
            champ_label.pack()
        bouton = tk.Button(self.f,text="OK", command=self.fetch)
        self.ent=tk.Entry(self.f)
        self.ent.pack()
        bouton.pack()
        self.f.mainloop()
        
        

        
    def mat_dotplot_selon_evalue(self):
        self.mat=np.zeros((len(self.cds1), len(self.cds2)))
        
        f=open(self.blast_file,'r')
        #mettre les accesions des cds homologues+evalue(du match) du fichier blast dans 
        #un dictionnaire de paires de gènes et leur e_value d'homologie
        
        #seuil=float(self.seuil)
        dic={}
        for line in f:
            if line[0]!='#':
                line=line.rstrip('/n')
                line=line.split()
                dic[(line[0], line[1])]=float(line[10])
            #mettre les accesions des cds des 2 génomes dans des listes(du file feature table)
            
        acc1=[]
        acc2=[]
    
        for line in self.cds1:
            acc1.append(line[10])
        for line in self.cds2:
            acc2.append(line[10])
                
   
        for pair in dic.keys():
        #print(dic[pair])
            if pair[0] in acc1 and pair[1] in acc2:
                #print(self.seuil)
                if dic[pair]<self.seuil:   #filtre: on ne prends que les matchs avec e_value inférieure au seuil
                    self.mat[acc1.index(pair[0])][acc2.index(pair[1])]=1
                             
                             
    def mat_dotplot_selon_id(self):
        self.mat=np.zeros((len(self.cds1), len(self.cds2)))
        f=open(self.blast_file,'r')
        
        #mettre les accesions des cds homologues+evalue(du match) du fichier blast dans 
        #un dictionnaire de paires de gènes et leur pourcentage d'identité %id
   
        dic={}
        for line in f:
            if line[0]!='#':
                line=line.rstrip('/n')
                line=line.split()
                dic[(line[0], line[1])]=float(line[2])
            #mettre les accesions des cds des 2 génomes dans des listes(du file feature table)
            
        acc1=[]
        acc2=[]    
        for line in self.cds1:
            acc1.append(line[10])
        for line in self.cds2:
            acc2.append(line[10])
                
   
        for pair in dic.keys():
         
            if pair[0] in acc1 and pair[1] in acc2:
                if dic[pair]>self.seuil:   #filtre: on ne prends que les matchs avec % d'id supérieure au seuil
                    self.mat[acc1.index(pair[0])][acc2.index(pair[1])]=1                      

                    
    def mat_dotplot_selon_couverture(self):
        self.mat=np.zeros((len(self.cds1), len(self.cds2)))
        f=open(self.blast_file,'r')
        dic={}   #dic pour couverture longueur de chaque gene_accession
        for line in f:
            if line[0]!='#':
                line=line.rstrip('/n')
                line=line.split()
  
                couverture1=float((float(line[7])-float(line[6]))/float(line[3]))
                couverture2=float((float(line[9])-float(line[8]))/float(line[3]))
                
                dic[(line[0], line[1])]=(couverture1, couverture2)
   
        acc1=[]
        acc2=[]    
        for line in self.cds1:
            acc1.append(line[10])
        for line in self.cds2:
            acc2.append(line[10])
         
        for pair in dic.keys():
            if pair[0] in acc1 and pair[1] in acc2:
                if float(((dic[pair][0]+dic[pair][1])/2)*100) > self.seuil:
                    self.mat[acc1.index(pair[0])][acc2.index(pair[1])]=1


    def dotplot(self):
        
        self.f2 = tk.Tk()
        
        self.f2.title("DotPlot")
        self.f2.geometry("500x500")
        if self.test=='e_value':
            self.mat_dotplot_selon_evalue()
        if self.test=="pourcentage d'id":
            self.mat_dotplot_selon_id()
        if self.test=="couverture":
            self.mat_dotplot_selon_couverture()
           
    
        print(self.seuil)
      
        ind=np.argwhere(self.mat==1)
        x=[]
        y=[]
        for i in  ind:
            x.append(i[0])
            y.append(i[1])
        self.x=x
        self.y=y

        fig = Figure(figsize = (5, 5), dpi = 100) 

        plot1 = fig.add_subplot(111) 

        plot1.scatter(x,y, c="black", s=2)
        canvas = FigureCanvasTkAgg(fig,master =self.f2)   
        canvas.draw()
        canvas.get_tk_widget().pack()
        self.f2.mainloop()
        
    def fetch(self):
        self.seuil=float(self.ent.get())
        self.f.destroy()
        self.dotplot()





class Premiere_Fenetre:
    def __init__(self,blast_file,cds1,cds2 ):
        
        self.blast_file=blast_file
        self.cds1= cds1
        self.cds2=cds2
        # Crée la fenêtre
        
        self.choix_test()
    def fonction_e_value(self):   
         
        self.fenetre_test.destroy()
        f=fenetres(self.blast_file,"e_value" ,self.cds1,self.cds2)
        
        
    def fonction_id(self):   
        
        self.fenetre_test.destroy()
        f=fenetres(self.blast_file,"pourcentage d'id" ,self.cds1,self.cds2)
        
    def fonction_couverture(self):   
       
        self.fenetre_test.destroy()
        f=fenetres(self.blast_file,"couverture" ,self.cds1,self.cds2)
      
       
    def choix_test(self):    
        self.fenetre_test = tk.Tk()

        self.fenetre_test.title("choisir le test")
        self.fenetre_test.geometry("300x200")
        boutG = tk.Button(self.fenetre_test,text="e_value", command=self.fonction_e_value)
        boutH = tk.Button(self.fenetre_test,text="pourcentage d'id", command=self.fonction_id)
        boutB = tk.Button(self.fenetre_test,text=" couverture", command=self.fonction_couverture)
        boutG.pack()
        boutH.pack()
        boutB.pack()
        self.fenetre_test.mainloop()
    

 
 
   
lancement=Premiere_Fenetre('QUERY-GCA_000009985.1_ASM998v1_protein__DB-GCA_000014865.1_ASM1486v1_protein.out', cds1, cds2)
  


