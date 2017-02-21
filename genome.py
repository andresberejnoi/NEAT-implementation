# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 11:54:58 2017

@author: andresberejnoi
"""
import Network
import random

#initialize the random seed
random.seed()
r = random.random

class GenomeError(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return msg
    __repr__ = __str__

class Genome(object):

    def __init__(self):
        self.connection_genes = {}

        #Node genes will be stored in separate dictionaries
        #in this way, selection for connection mutations will be easier
        self.node_genes = {}                        #a dictionary of the hidden nodes in the network
        self._input_nodes = {}
        self._output_nodes = {}
        self._hidden_nodes = {}

        self.current_connections = set()
        self.max_nodeID = 0
        self.fitness = 0
        #self.connection_ids = list(self.connection_genes)
        #self._update_id_lists()

    def _update_sets(self):
        #this could be done in a more efficient way
        for ID in self.connection_ids:
            g =self.connection_genes[ID]
            self.current_connections.add(frozenset((g.sendingNode,g.receivingNode)))        #creates a set of sets, and therefore no pair of nodes can be repeated
    def _update_id_lists(self):
        """Updates the lists that keep track of the different IDs"""
        self.connection_ids = list(self.connection_genes)
        self.node_ids = list(self.node_genes)
        self.max_nodeID = max(self.node_genes)

    #operator overloading
    '''
    def __getitem__(self, key):
        """Allows for dictionary or indexing behaviour"""
        return self.genes[key]
    def __setitem__(self,key,val):
        self.genes[key] = val
    '''

    def add_connection_gene(self, genes):
        """TODO: allow this function to take a list of genes as argument, so it can be called fewer times"""
        #First, if the argument is a single gene check if it is a Genome object
        if isinstance(genes, Genome):
            ID = gene.get_id()      #for a connection gene, the ID will be the innovation number
            try:
                if id in self.connection_genes:
                    raise GenomeError
                else:
                    self.connection_genes[ID] = gene
                    self.current_connections.add(frozenset((gene.sendingNode,gene.receivingNode)))

            except GenemeError:
                print(("""Gene with id {0} already in genome""".format(ID)))
        elif isinstance(genes, (list,tuple,dict)):
            for gene in genes:
                ID = gene.get_id()
                try:
                    if ID in self.connection_genes:
                        raise GenomeError
                    else:
                        self.connection_genes[ID] = gene
                        self.current_connections.add(frozenset((gene.sendingNode,gene.receivingNode)))
                except GenemeError:
                    print(("""Gene with id {0} already in genome""".format(ID)))


        self._update_id_lists()

    def add_node_gene(self, genes):
        """TODO: allow this function to take a list of genes as argument, so it can be called fewer times"""
        if isinstance(genes, NodeGene):
            ID = gene.get_id()
            try:
                if id in self.node_genes:
                    raise GenomeError
                else:
                    self.node_genes[ID] = gene
            except GenemeError:
                print(("""Gene with id {0} already in genome""".format(ID)))

        elif isinstance(genes, (list, tuple,dict)):
            for gene in genes:
                ID = gene.get_id()
                if ID in self.node_genes:
                    raise GenomeError("""ID {0} already exists in genome""".format(ID))
                else:
                    self.node_genes[ID] = gene
        self._update_id_lists()

    def mutate(self, mRate=0.01):
        """
        TODO: this function should return a new genome without modifying the current one

        mRate: the mutation rate. This should probably come from a configuration file later on
        There are 3 types of possible mutations:
        Weight mutation: every weight in the network can change by some amount
        Mutate add connection: two previously unconnected nodes are connected by a new link.
                                This adds a connection gene to the genome
        Mutate add node: an existing connection is split and a node is created in the middle.
                        This adds 2 new connection genes and one node gene to the genome
        """
        add_connection_mRate = 0.2
        add_node_mRate = 0.2
        change_weight_mRate = 0.2
        r = random.random           #function to generate a random number between 0 and 1
        if r()<add_node_mRate:           #there should actually be different mutation rates for the 3 different kinds of mutations possible
            self._mutate_add_node()
        elif r()<add_connection_mRate:
            self._mutate_add_connection()
        else:
            #This mutation modifies an existing connection
            self._mutate_modify_connection()

        #We should also keep track of the structural innovations in a global list so that if two mutations result in the same innovation, they receive the same innovation number

        self._update_id_lists()

    def _mutate_add_connection(self):
        """A connection will be added between any two unconnected nodes"""
        #determine which nodes are not currently connected
        #select two node ids at random and check that they do not have a connection
        n1 = random.choice(self.node_ids)           #modify this so that the input node cannot be chosen from the output nodes of the network
        n2 = random.choice(self.node_ids)           #modify this so that the output node cannot be chosen from the input nodes of the network

        while set(n1,n2) in self.current_connections:
            n1 = random.choice(self.node_ids)
            n2 = random.choice(self.node_ids)

        #create a new connection gene
        new_gene = WeightGene(sendingNode=n1,receivingNode=n2,weight=1.0,enabled=1)
        self.current_connections.add(frozenset(n1,n2))

    def _mutate_add_node(self):
        """A random connection will be split into two and a new node will be created in between"""

        #first, create a new node gene:
        ID = self.max_nodeID + 1                       #The new node will have an ID higher than the previous higher
        new_node = NodeGene(ID=ID, nodeType="H", actFun="tanh")

        #then, select an existing connection at random to split:
        c = random.choice(self.connection_ids)            #I read somewhere that random.choice is slow, so maybe I should look at another option
        new_c1, new_c2 = c.split(middleNode=ID)      #these two are the new connections

        #Add all three new genes to the genome
        self.add_connection_gene(new_c1)
        self.add_connection_gene(new_c2)
        self.add_node_gene(new_node)

    def _mutate_modify_connection(self):
        #select a random weight to modify
        #A mutation can be a change in weight, toggling the enable bit
        if r() < 0.8:           # there is an 80% chance that the mutation will change the weight (80% is just arbitrary)
            pass


    def crossover(self, genome2):
        """
        Performs crossover between self and genome2 and returns a new genome. The code below can be simplified a lot more.
        """
        new_genome = Genome()
        match_genes = []
        non_match_genes_1 = []
        non_match_genes_2 = []

        IDs = zip(self.connection_ids, genome2.connection_ids)
        for id1,id2 in IDs:
            if id1 in genome2.connection_ids:
                match_genes.append(id1)
            else:
                non_match_genes_1.append(id1)
                non_match_genes_2.append(id2)

        #Crossover matching connection genes
        for ID in match_genes:
            # For matching genes, the new gene will be copied as is, from one of the parents at random:
            if r() < 0.5:
                new_gene = self.connection_genes[ID].copy()
            else:
                new_gene = genome2.connection_genes[ID].copy()
            new_genome.add_connection_gene(new_gene)

            '''
            g1 = self.connection_genes[ID]
            g2 = genome2.connection_genes[ID]
            new_gene = g1.reproduce(g2)
            new_genome.add_connection_gene(new_gene)         # this needs to change to actually adding a gene to a genome object instead of a list
            '''

        #Now deal with the non-matching connection genes, and with node genes
        f1 = self.fitness
        f2 = genome2.fitness

        if f1 > f2:
            #Adding connection genes
            for ID in non_match_genes_1:
                new_gene = self.connection_genes[ID].copy()
                new_genome.add_connection_gene(new_gene)
        elif f1 < f2:
            for ID in non_match_genes_2:
                new_gene = genome2.connection_genes[ID].copy()
                new_genome.add_connection_gene(new_gene)
        else:
            for ID in non_match_genes_1:
                new_genome.add_connection_gene(self.connection_genes[ID].copy())
            for ID in non_match_genes_2:
                new_genome.add_connection_gene(genome2.connection_genes[ID].copy())

        #Now perform crossover of node genes:
        #Adding node genes
        for nodeID in self.node_ids:
            if nodeID in genome2.node_ids:
                new_node_gene = self.node_genes[nodeID].reproduce(genome2.node_genes[nodeID])
                new_genome.add_node_gene(new_node)
            else:
                #here we deal with non-matching nodes
                if f1>f2:
                    for nodeID in self.node_ids:
                        new_node_gene = self.node_genes[nodeID].copy()
                        new_genome.add_node_gene(new_node_gene)
                elif f2>f1:
                    for nodeID in genome2.node_ids:
                        new_node_gene = genome2.node_genes[nodeID].copy()
                        new_genome.add_node_gene(new_node_gene)
                else:
                    #if both genomes have the same fitness, then all the nodes are inherited.
                    for nodeID in self.node_ids:
                        new_node_gene = self.node_genes[nodeID].copy()
                        new_genome.add_node_gene(new_node_gene)
                    for nodeID in genome2.node_ids:
                        new_node_gene = genome2.node_genes[nodeID].copy()
                        new_genome.add_node_gene(new_node_gene)

        return new_genome



class WeightGene(object):
    def __init__(self, sendingNode, receivingNode, weight, enabled, hisMark=None):
        """
        """
        self.sendingNode = sendingNode
        self.receivingNode = receivingNode
        self.weight = weight
        self.enabled = enabled
        #Update the global innovation number
        if hisMark==None:
            global INNOVATION_NUMBER
            INNOVATION_NUMBER += 1
            self.hisMark = INNOVATION_NUMBER           # the historical marking (innovation number; the one mentioned in the NEAT papers)
        else: self.hisMark = hisMark

    #Getters:
    def get_id(self):
        return self.hisMark

    def checkHomology(self, gene2):
        """
        Checks for homology by comparing the innovation number in gene2 with self.
        Returns True if the two genes are related and False if they are not
        """
        if self.hisMark == gene2.hisMark:
            return True
        return False

    def copy(self):
        """Creates an identical copy of the current gene"""
        geneCopy = WeightGene(self.sendingNode, self.receivingNode, self.weight, self.enabled, self.hisMark)
        return geneCopy

    def split(self, middleNode):
        """
        middleNode: int; this is the node id that will be in the middle of the two new connections
        Split the connection into 2 different ones:
            The 1st connection begins at self.sendingNode and ends at middleNode; The weight will be 1.0
            The 2nd connection begins at self.middleNode and ends at self.receivingNdoe; the weight will be equal to the original connection weight
        The current connection will be disabled.
        """
        connection1 = WeightGene(sendingNode=self.sendingNode, receivingNode=middleNode, weight=1.0, enabled=1)           # how should I increase the historical marking?
        connection2 = WeightGene(sendingNode=middleNode,receivingNode=self.receivingNode, weight=self.weight, enabled=1)        # Also here, the histocal marking needs to be increased

        self._disable_gene()
        return connection1, connection2
    def reproduce(self, gene2):
        """
        Creates a new gene that inherits attributes from each parent gene.
        """
        #TODO: automate the process of selecting attributes (look into dir() and getattr()
        #r = random.random

        #TODO: Since genes with the same innovation number will have the same input/output node,
        #   then it is not necessary to randomly choose between the two
        #sendNode_inherit = self.sendingNode if r() < 0.5 else gene2.sendingNode
        #eceiveNode_inherit = self.receivingNode if r() < 0.5 else gene2.receivingNode
        sendNode_inherit = self.sendingNode
        receiveNode_inherit = self.receivingNode

        weight_inherit = self.weight if r() < 0.5 else gene2.weight
        hisMark_inherit = self.hisMark                  # this mark should be the same in both genes (I think), so there is no need to select randomly
        #To decide on the enabled bit, if both parent genes are enabled, then the offspring will be enabled.
        #   if one of them is disabled, then there will be a random chance that the offspring is disabled
        #   if both parent genes are disabled, the offspring will be disabled (I am not sure if this part follows the NEAT algorithm)
        if self.enabled + gene2.enabled == 2:
            enabled_inherit = 1
        elif self.enabled + gene2.enabled == 1:
            enabled_inherit = self.enabled if r() < 0.5 else gene2.enabled
        else:
            enabled_inherit = 0

        newGene = WeightGene(sendNode_inherit, receiveNode_inherit, weight_inherit, enabled_inherit, hisMark_inherit)
        return newGene

    def _disable_gene(self):
        self.enabled = 0

    def _toggle_enable(self):
        """Toggles the enable bit"""
        self.enabled = abs(self.enabled-1)

    def _mutate_weight(self):
        """"""
    def _mutate_connection(self):
        """"""
    def _mutate_enablebit(self):
        """"""
    def mutate(self):
        """
        Brings in all the mutation functions in one place
        """


class NodeGene(object):
    def __init__(self, ID, nodeType, actFun):
        """
        ID: int; the id associated with the node
        nodeType: string; string indicating the node type ('I' for input, 'H' for hidden, 'O' for output)
        actFun: string; the name of the activation function for the neuron. It will be looked up in a funciton dictionary (e.g. 'tanh', 'sigmoid')
        """
        self.id = ID
        #self.inputs = inputs        #a list of the neurons that are inputs to this one
        #self.outputs = outputs      #a list of the neurons that are outputs to this one
        self.nodeType = nodeType
        self.actFun = actFun        #the activation function for this neuron

    #Getters:
    def get_id(self):
        return self.id

    def copy(self):
        """
        Create a new NodeGene that has all the same properties as the reproducing gene
        """
        geneCopy = NodeGene(self.id, type, self.actFun, )
        return geneCopy

    def reproduce(self, geneParent2):
        """
        Create a new NodeGene that is a combination of this gene and geneParent2
        Traits from each parent gene will be passed on byy chance
        """
        #TODO: use an automated way of looping through the attributes instead of handwriting this part (the combination of dir() and getattr()
        r = random.random
        id_inherit = self.id if r() < 0.5 else geneParent2.id           #I think this is not the correct way to do this
        type_inherit = self.nodeType if r() < 0.5 else geneParent2.nodeType
        actFun_inherit = self.actFun if r() < 0.5 else geneParent2.actFun

        newGene = NodeGene(id_inherit, type_inherit, actFun_inherit)

        return newGene


