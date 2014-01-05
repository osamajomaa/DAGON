"""
a simple Gene Ontology represented by a Directed Asyclic Graph with 
capabilities like adding terms, relations, inferencing and drawing the ontology"
"""

import networkx as nx # NetworkX: the main library used to create and manipulate the DAG
import matplotlib.pyplot as plt # matplotlib: the library used to draw GO

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
        self.add_term("GO:0000001", "cellular_component", Aspect.cc, "cellular component aspect of the GO")
        self.add_term("GO:0000002", "biological_process", Aspect.bp, "biological process aspect of the GO")
        self.add_term("GO:0000003", "molecular function", Aspect.mf, "molecular function aspect of the GO")
        
    def add_term(self, uid, tname, tnamespace, tdefinition):
        """
        Again, this is a simple GO, we added only the essential attributes of a term; the unique id, name, namespace or aspect and the definition.
        Real GO has other optional attributes.
        """
        
        self.Ontology.add_node(uid, name=tname, namespace=tnamespace, definition=tdefinition)    
        
    def add_relation(self, term1, term2, rel):
        
        # when creating a relation between two terms, first of all they should be in the same subontology(aspect)
        if self.Ontology.node[term1]['namespace'] != self.Ontology.node[term2]['namespace']:
            print("The relation was not added as the two terms belong to different aspects of GO")
            return    
        self.Ontology.add_edge(term1, term2, relation=rel)
        #The new relation must not form a cycle in the directed graph
        if (nx.algorithms.dag.is_directed_acyclic_graph(self.Ontology) == False):
            self.Ontology.remove_edge(term1, term2)
            print("Relation was not added as it would create cycles in the graph!")
        else:
            print("Relation was added successfully!")
    
    def Infer_sem_related_terms(self,term):
        """
        Inferencing over the ontology is one of the key features of GO. 
        Starting from a term we can get all terms which are semantically related to this term by either is_a or part_of relation.
        """
        
        related_terms = {}
        self.is_a_inference(term,related_terms)
        is_a_terms = list(related_terms.keys())
        
        related_terms = {}
        self.part_of_inference(term,related_terms)
        part_of_terms = list(related_terms.keys())
        
        return is_a_terms, part_of_terms
        
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
    
    def draw_GO(self):
        """
        Draw GO as a DAG
        """
        nx.draw(self.Ontology)
        plt.show()


#Sample code to test GO class

go = GO()
go.add_term("GO:0000004", "1", Aspect.mf, "")
go.add_relation("GO:0000004", "GO:0000003", Relation.is_a)

go.add_term("GO:0000005", "2", Aspect.mf, "")
go.add_relation("GO:0000005", "GO:0000003", Relation.is_a)

go.add_term("GO:0000006", "0", Aspect.mf, "")
go.add_relation("GO:0000004", "GO:0000006", Relation.part_of)

go.add_term("GO:0000007", "3", Aspect.mf, "")
go.add_relation("GO:0000007", "GO:0000004", Relation.is_a)

go.add_term("GO:0000008", "4", Aspect.mf, "")
go.add_relation("GO:0000008", "GO:0000004", Relation.is_a)
go.add_relation("GO:0000008", "GO:0000005", Relation.part_of)

go.add_term("GO:0000009", "5", Aspect.mf, "")
go.add_relation("GO:0000009", "GO:0000005", Relation.part_of)

go.add_term("GO:0000010", "6", Aspect.mf, "")
go.add_relation("GO:0000010", "GO:0000007", Relation.is_a)

go.add_term("GO:0000011", "7", Aspect.mf, "")
go.add_relation("GO:0000011", "GO:0000004", Relation.is_a)
go.add_relation("GO:0000011", "GO:0000008", Relation.part_of)

go.add_term("GO:0000012", "8", Aspect.mf, "")
go.add_relation("GO:0000012", "GO:0000009", Relation.is_a)

is_a_terms, part_of_terms = go.Infer_sem_related_terms("GO:0000007")
print("Set of GO terms resulted by inferencing over GO which are semantically related to GO:0000007 by is_a relation")
print (is_a_terms)
print("======================")
print("Set of GO terms resulted by inferencing over GO which are semantically related to GO:0000007 by part_of relation")
print(part_of_terms)

go.draw_GO()
