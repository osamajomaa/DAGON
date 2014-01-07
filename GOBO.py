"""
a simple Gene Ontology represented by a Directed Asyclic Graph with 
capabilities like adding terms, relations, inferencing and drawing the ontology"
"""

import networkx as nx # NetworkX: the main library used to create and manipulate the DAG
import go_obo_parser as gop
import matplotlib.pyplot as plt # matplotlib: the library used to draw GO
import pickle

def enum(**enums): # since python 3.3 does not have an enum type, I created one 
    return type('Enum', (), enums)

Relation = enum(is_a="is_a", part_of="part_of") # Relation types in our simple GO. Real GO has more relation types like regulates
Aspect = enum(cc="cellular_component", bp="biological_process", mf="molecular_function") # The three aspects of GO  


class GO:    
    
    def __init__(self):
        """
        We start off by creating the three basic terms which all other terms connect to forming the three aspects of go.
        These three terms are cellular component, biological process and molecular function.        
        """        
        self.Ontology = nx.classes.DiGraph() 
        
    def __get_aspect__(self, aspect):
        if aspect == "cellular_component":
            return Aspect.cc
        elif aspect == "biological_process":
            return Aspect.bp
        else:
            return Aspect.mf
        
    def upload_go_from_dagon(self, dagon_file):
        d_file = open(dagon_file)
        self.Ontology = pickle.load(d_file)        
    
    def upload_go_from_obo(self, go_obo_file):
        """"
        Upload an OBO file into the GO data structure. 
        We upload all terms then all relations because most relations have terms that hasn't been 
        added to the GO structre yet.
        """
        # Firest, upload all terms in the obo file 
        for go_term in gop.parseGOOBO(go_obo_file):
            if not go_term.has_key('is_obsolete'):
                print go_term["id"][0]
                go.add_term(go_term["id"][0], go_term["name"][0], self.__get_aspect__(go_term["namespace"][0]), go_term["def"][0])  
        
        # Second, upload all relations in the obo file        
        for go_term in gop.parseGOOBO(go_obo_file):
            print i
            if go_term.has_key('is_a'):
                for parent in go_term['is_a']:
                    self.add_relation(go_term["id"][0], parent["id"], Relation.is_a)
            if go_term.has_key("relationship"):
                for parent in go_term['relationship']:
                    if parent[:7] == "part_of":
                        self.add_relation(go_term["id"][0], parent[8:18], Relation.part_of)
            
            # serialize the ontology for a faster retrieval later
            fHandler = open("go.dagon", 'wb')
            pickle.dump(go.Ontology, fHandler)
            fHandler.close()            
        
    def add_term(self, uid, tname, tnamespace, tdefinition):
        """
        Again, this is a simple GO, we added only the essential attributes of a term; the unique id, name, namespace or aspect and the definition.
        Real GO has other optional attributes.
        """        
        self.Ontology.add_node(uid, name=tname, namespace=tnamespace, definition=tdefinition)    
        
    def add_relation(self, term1, term2, rel): 
        self.Ontology.add_edge(term1, term2, relation=rel)
    
    def get_ancestors(self,term, rel):
        """
        Get all ancestors for a term that are semantically related through part_of or is_a relation        
        """        
        ancestors = {}
        if rel == Relation.is_a:
            self.is_a_inference(term,ancestors)
        elif rel == Relation.part_of:
            self.part_of_inference(term,ancestors)
        
        return list(ancestors.keys())
        
    def is_a_inference(self, term, is_a_terms):
        """
        A is_a B and B is_a C --> A is_a C
        """
        for nbr in self.Ontology.neighbors_iter(term):
            if self.Ontology[term][nbr]['relation'] == Relation.is_a: 
                is_a_terms[nbr] = None
                self.is_a_inference(nbr,is_a_terms)
    
    def part_of_inference(self, term, part_of_terms):
        """
        A is_a B and B part_of C --> A part_of C
        A part_of B and B is_a C --> A part_of C
        """
        for nbr in self.Ontology.neighbors_iter(term):
                if self.Ontology[term][nbr]['relation'] == Relation.part_of:
                    part_of_terms[nbr] = None
                    self.part_of_inference_Helper(nbr, part_of_terms)
                elif self.Ontology[term][nbr]['relation'] == Relation.is_a:
                    self.part_of_inference(nbr,part_of_terms)    
    
    def part_of_inference_Helper(self, term, part_of_terms):
        for nbr in self.Ontology.neighbors_iter(term):
                if self.Ontology[term][nbr]['relation'] == Relation.is_a or self.Ontology[term][nbr]['relation'] == Relation.part_of: 
                    part_of_terms[nbr] = None
                    self.part_of_inference_Helper(nbr,part_of_terms)
    
    def save_ancestors(self, term, ancestors, rel):
        """"
        Save ancestors in a conveinent format in txt file
        """
        fname = rel+"_ancestors_file.txt"
        anc_file = open(fname, "w")
        txt = "This file contains all the ancestors of the GO term:\n"
        txt += "id: "+ term + "\n"
        txt += "name: " + self.Ontology.node[term]["name"] + "\n"
        txt += "namespace: " + self.Ontology.node[term]["namespace"] + "\n"
        txt += "def: " + self.Ontology.node[term]["definition"]  + "\n"
        txt += "that are related by " + rel + " relationship.\n\n"
        txt += "Ancestors:\n"
        for anc in ancestors:
            txt += "id: "+ anc + "\n"
            txt += "name: " + self.Ontology.node[anc]["name"] + "\n"
            txt += "namespace: " + self.Ontology.node[anc]["namespace"] + "\n"
            txt += "def: " + self.Ontology.node[anc]["definition"]  + "\n\n"            
        anc_file.write(txt)
    
    def draw_GO(self):
        """
        Draw GO as a DAG
        """
        nx.draw(self.Ontology)
        plt.show()


#Sample code to test GO class

go = GO()
#go.upload_go_from_obo("go.obo")
go.upload_go_from_dagon("go.dagon")
part_of_ancestors = go.get_ancestors("GO:0000083", Relation.part_of)
is_a_ancestors = go.get_ancestors("GO:0000056", Relation.is_a)
go.save_ancestors("GO:0000083", part_of_ancestors,Relation.part_of)
go.save_ancestors("GO:0000056", is_a_ancestors,Relation.is_a)