import json
import csv
import numpy as np

class Dat:
    
    def __init__(self):
        self.pathdict = {
            'wiki' : '/Users/kylesargent/Desktop/GLoVe + VerbNet/glove.6B.300d.txt', 
            'twitter' : '/Users/kylesargent/Desktop/GLoVe + VerbNet/glove.twitter.27B.100d.txt'
        }
        
        self.gen_dict = {}
        self.GLoVe_dict_list = []
        self.frame_list = frame_list = ['It V', 'It V NP.theme', 'It V PP.experiencer that S', 'It V PP.theme', 'It V that S', 'NP NP V ADVP-Middle together', 'NP NP V together', 'NP v', 'NP V', 'NP V', 'NP V ADJ', 'NP V ADJ PP.experiencer', 'NP V ADJ-Middle', 'NP V ADJP NP', 'NP V ADV', 'NP V ADV-Middle', 'NP V ADV-Middle PP', 'NP V ADV-Middle PP.location', 'NP V ADV-Middle together', 'NP V ADVP', 'NP V ADVP-Middle', 'NP V ADVP-Middle PP', 'NP V apart', 'NP V down NP', 'NP V for NP S_INF', 'NP V how S', 'NP V how S_INF', 'NP V NP', 'NP V NP', 'NP V NP ADJ', 'NP V NP ADJ PP.instrument', 'NP V NP ADJP', 'NP V NP ADJP PP.instrument', 'NP V NP ADJP PP.result', 'NP V NP ADJP-Result', 'NP V NP ADV', 'NP V NP ADVP', 'NP V NP apart', 'NP V NP down', 'NP V NP how S', 'NP V NP how S_INF', 'NP V NP NP', 'NP V NP NP PP.asset', 'NP V NP NP PP.theme', 'NP V NP NP together', 'NP V NP NP.asset', 'NP V NP P.asset', 'NP V NP PP', 'NP V NP PP PP.instrument', 'NP V NP PP.asset', 'NP V NP PP.attribute', 'NP V NP PP.attribute PP.extent', 'NP V NP PP.beneficiary', 'NP V NP PP.co-agent', 'NP V NP PP.co-patient', 'NP V NP PP.co-theme', 'NP V NP PP.destination', 'NP V NP PP.destination PP.initial_location', 'NP V NP PP.destination-Conative', 'NP V NP PP.destinations', 'NP V NP PP.goal', 'NP V NP PP.initial_location', 'NP V NP PP.initial_location PP.destination', 'NP V NP PP.instrument', 'NP V NP PP.location', 'NP V NP PP.material', 'NP V NP PP.material PP.asset', 'NP V NP PP.material PP.beneficiary', 'NP V NP PP.patient', 'NP V NP PP.predicate', 'NP V NP PP.product PP.beneficiary', 'NP V NP PP.recipient', 'NP V NP PP.recipient PP.asset', 'NP V NP PP.recipient PP.theme', 'NP V NP PP.result', 'NP V NP PP.result PP.instrument', 'NP V NP PP.source', 'NP V NP PP.source NP.asset', 'NP V NP PP.source PP.beneficiary', 'NP V NP PP.source PP.instrument', 'NP V NP PP.source S', 'NP V NP PP.source S_ING', 'NP V NP PP.source whether S', 'NP V NP PP.stimulus', 'NP V NP PP.theme', 'NP V NP PP.theme', 'NP V NP PP.topic', 'NP V NP PP.topic what S', 'NP V NP PP.value', 'NP V NP S', 'NP V NP S-INF', 'NP V NP S-Quote', 'NP V NP S_INF', 'NP V NP S_ING', 'NP V NP that S', 'NP V NP to be ADJ', 'NP V NP to be NP', 'NP V NP together', 'NP V NP up', 'NP V NP what S', 'NP V NP what S_INF', 'NP V NP whether S', 'NP V NP whether S_INF', 'NP V NP-ATTR-POS', 'NP V NP-dative NP', 'NP V NP-Dative NP', 'NP V NP-Fulfilling PP', 'NP V NP-Fulfilling PP.theme', 'NP V NP-PRO-ARB', 'NP V NP.attribute', 'NP V NP.beneficiary NP', 'NP V NP.beneficiary NP PP', 'NP V NP.destination', 'NP V NP.destination PP.theme', 'NP V NP.experiencer', 'NP V NP.initial_location', 'NP V NP.location', 'NP V NP.location PP.theme', 'NP V NP.material', 'NP V NP.material PP.product', 'NP V NP.patient', 'NP V NP.patient PP.material PP.result', 'NP V NP.patient PP.result', 'NP V NP.product PP.material', 'NP V NP.recipient', 'NP V NP.recipient how S_INF', 'NP V NP.recipient S', 'NP V NP.recipient S-Quote', 'NP V NP.recipient S_INF', 'NP V NP.recipient that S', 'NP V NP.recipient when S_INF', 'NP V NP.source', 'NP V NP.stimulus', 'NP V NP.theme', 'NP V NP.theme (PP)', 'NP V NP.theme PP', 'NP V NP.theme PP.destination', 'NP V NP.theme PP.instrument', 'NP V NP.topic', 'NP V NP.value', 'NP V NP:', 'NP V out', 'NP V PP', 'NP V PP ADV-Middle', 'NP V PP PP', 'NP V PP-Conative', 'NP V PP.attribute', 'NP V PP.attribute S_ING', 'NP V PP.beneficiary', 'NP V PP.co-agent', 'NP V PP.co-agent how S', 'NP V PP.co-agent how S_INF', 'NP V PP.co-agent PP.goal', 'NP V PP.co-agent PP.theme', 'NP V PP.co-agent PP.theme S_ING', 'NP V PP.co-agent PP.theme what S', 'NP V PP.co-agent PP.theme what S_INF', 'NP V PP.co-agent PP.theme whether S', 'NP V PP.co-agent PP.theme whether S_INF', 'NP V PP.co-agent PP.topic', 'NP V PP.co-agent PP.topic what S', 'NP V PP.co-agent PP.topic what S_INF', 'NP V PP.co-agent PP.topic whether S_INF', 'NP V PP.co-agent PP.topic whether/if S', 'NP V PP.co-patient', 'NP V PP.co-theme', 'NP V PP.destination', 'NP V PP.destination NP', 'NP V PP.experiencer', 'NP V PP.extent', 'NP V PP.goal', 'NP V PP.goal what S_INF', 'NP V PP.initial_loc', 'NP V PP.initial_loc PP.destination', 'NP V PP.initial_location', 'NP V PP.initial_location PP.destination', 'NP V PP.instrument', 'NP V PP.location', 'NP V PP.location PP.theme', 'NP V PP.material PP.result', 'NP V PP.patient', 'NP V PP.recipient', 'NP V PP.recipient how S', 'NP V PP.recipient how S_INF', 'NP V PP.recipient NP', 'NP V PP.recipient PP.topic', 'NP V PP.recipient PP.topic what S', 'NP V PP.recipient PP.topic what S_INF', 'NP V PP.recipient PP.topic whether S_INF', 'NP V PP.recipient S-Quote', 'NP V PP.recipient S_INF', 'NP V PP.recipient that S', 'NP V PP.recipient what S', 'NP V PP.recipient what S_INF', 'NP V PP.recipient when S_INF', 'NP V PP.recipient whether S', 'NP V PP.recipient whether/if S', 'NP V PP.result', 'NP V PP.result NP', 'NP V PP.result PP.instrument', 'NP V PP.source', 'NP V PP.source how S', 'NP V PP.source how S_INF', 'NP V PP.source PP.goal', 'NP V PP.source PP.theme what S', 'NP V PP.source PP.theme what S_INF', 'NP V PP.source PP.theme whether S_INF', 'NP V PP.source PP.theme whether/if S', 'NP V PP.source S_ING', 'NP V PP.stimulus', 'NP V PP.stimulus how/whether S', 'NP V PP.stimulus S_ING', 'NP V PP.stimulus what S', 'NP V PP.theme', 'NP V PP.theme NP S_ING', 'NP V PP.theme NP.location', 'NP V PP.theme PP.location', 'NP V PP.theme PP.source', 'NP V PP.theme S', 'NP V PP.theme S_INF', 'NP V PP.theme S_ING', 'NP V PP.theme what S', 'NP V PP.theme what S_INF', 'NP V PP.theme whether S_INF', 'NP V PP.theme whether/if S', 'NP V PP.time', 'NP V PP.topic', 'NP V PP.topic NP S_ING', 'NP V PP.topic PP.co-agent', 'NP V PP.topic PP.recipient', 'NP V PP.topic S_ING', 'NP V PP.topic what S', 'NP V PP.topic what S_INF', 'NP V PP.topic whether S', 'NP V PP.topic whether S_INF', 'NP V PP.trajectory', 'NP V PP.value', 'NP V S', 'NP V S-Quote', 'NP V S_INF', 'NP V S_ING', 'NP V that S', 'NP V that S PP.theme', 'NP V together', 'NP V together ADV-Middle', 'NP V up NP', 'NP V what S', 'NP V what S_INF', 'NP V when S_INF', 'NP V whether S', 'NP V whether S_INF', 'NP V whether/if S', 'NP V whether/if S_INF', 'NP V why S', 'NP.agent V', 'NP.agent V NP', 'NP.asset V NP', 'NP.asset V NP NP', 'NP.asset V NP PP', 'NP.asset V PP', 'NP.attribute V', 'NP.attribute V NP.extent', 'NP.attribute V PP.extent', 'NP.cause V NP', 'NP.instrument V ADVP', 'NP.instrument V NP', 'NP.location V', 'NP.location V NP', 'NP.location V NP.theme', 'NP.location V PP.agent', 'NP.location V PP.theme', 'NP.material V NP', 'NP.material V PP.product', 'NP.patient V', 'NP.patient V PP.attribute', 'NP.patient V PP.material PP.result', 'NP.patient V PP.result', 'NP.product V PP.material', 'NP.theme V', 'NP.theme V NP', 'NP.theme V PP.location', 'NP.theme V PP.source', 'Passive', 'PP.location there V NP', 'PP.location V NP', 'PP.location V PP.theme', 'That S.stimulus V', 'There V NP', 'There V NP PP', 'There V NP PP.location', 'There V PP NP' ]

        self.verbnet_dict = {}
        reader = csv.reader(open('/Users/kylesargent/Desktop/GLoVe + VerbNet/verb_frames.csv', 'r'))
        for row in reader:
            k, v = row
            self.verbnet_dict[k] = json.loads(v)
            
    def __write_dict(self, path):
        dct = {}
        file = open(path)
        for row in file:
            key = row[:row.index(' ')]
            value = list(map(lambda x : float(x), str.split(row[row.index(' ') + 1:])))
            dct[key] = value
        return(dct)
            
    def add_dict(self, pathname):
        self.GLoVe_dict_list += [self.__write_dict(self.pathdict[pathname])]

    def make_gen_dict(self):
        def safe_index(dct, key): 
            try: return(dct[key]) 
            except KeyError: return(False)
        self.gen_dict = {}
        for verb, verbnet_vec in self.verbnet_dict.items():
            l = list(map(lambda dct : safe_index(dct, verb), self.GLoVe_dict_list))
            if not False in l:
                self.gen_dict[verb] = (verbnet_vec, np.concatenate(l))
        
    def getXy(self):
        dat = list(self.gen_dict.items())
        X = np.array([value[1] for _key, value in dat])
        y = np.array([value[0] for _key, value in dat])
        return(X,y)                         