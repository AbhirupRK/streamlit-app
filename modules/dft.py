import numpy as np

class ReadVaspOutput:
    def __init__(self, file_content):
        # If file_content is a path (for local testing), open it
        if isinstance(file_content, str):
            with open(file_content, "r") as file:
                self.lines = file.readlines()
        # If file_content is a file-like object (from Streamlit), read it
        elif hasattr(file_content, 'read'): # Check if it has a read method
            # Read as text (assuming UTF-8 encoding, adjust if needed)
            self.lines = file_content.read().decode("utf-8").splitlines()
        else:
            raise TypeError("file_content must be a file path or a file-like object.")

        self.data = {}
        self._parse_file()
        self._validate_keys()
    
    def _parse_file(self):
        lines = self.lines
        species = []

        for line in lines:
            if "VRHFIN" in line:
                species.append(line.split("=")[1].split(":")[0].strip())
            
        for i in range(len(lines)):
            if "LSORBIT" in lines[i]:
                self.data["soc"] = list(map(str, lines[i].split()))[2]
            if "ions per type =" in lines[i]:
                nions = [int(x) for x in list(map(str,lines[i].split()))[4:]]
                self.data["nions"] = nions; natoms = sum(nions)
            if "direct lattice vectors" in lines[i]:
                lattice = []
                for j in range(i+1, i+4):
                    lattice.append(list(map(float, lines[j].split()))[:3])
                self.data["lattice"] = np.array(lattice)
            if "position of ions in cartesian coordinates  (Angst):" in lines[i]:
                pos_cart = []
                for j in range(i+1, i+1+natoms):
                    pos_cart.append(list(map(float, lines[j].split()))[:3])
                self.data["pos_cart"] = np.array(pos_cart)
            if "magnetization (x)" in lines[i]:
                m_x = []
                for j in range(i+4,i+4+natoms):
                    m_x.append(list(map(float, lines[j].split()))[1:])
                m_x = np.array(m_x)
            if "magnetization (y)" in lines[i]:
                m_y = []
                for j in range(i+4,i+4+natoms):
                    m_y.append(list(map(float, lines[j].split()))[1:])
                m_y = np.array(m_y)
            if "magnetization (z)" in lines[i]:
                m_z = []
                for j in range(i+4,i+4+natoms):
                    m_z.append(list(map(float, lines[j].split()))[1:])
                m_z = np.array(m_z)

        self.data["species"] = species
        if self.data["soc"] == "F":
            m_z = m_x; m_x = np.zeros((natoms,4)); m_y = np.zeros((natoms,4))
        self.data["magmom"] = np.stack((m_x[:,3],m_y[:,3],m_z[:,3]), axis=1)
    
    def _validate_keys(self):
        required_keys = ["lattice", "pos_cart", "nions", "species", "magmom"]
        missing_keys = [key for key in required_keys if key not in self.data]
        if missing_keys:
            raise ValueError(f"Missing required keys in data: {missing_keys}")

    def get(self, key, default=None):
        return self.data.get(key, default)
    

class ReadSiestaOutput:
    def __init__(self, file_content):
        # If file_content is a path (for local testing), open it
        if isinstance(file_content, str):
            with open(file_content, "r") as file:
                self.lines = file.readlines()
        # If file_content is a file-like object (from Streamlit), read it
        elif hasattr(file_content, 'read'): # Check if it has a read method
            # Read as text (assuming UTF-8 encoding, adjust if needed)
            self.lines = file_content.read().decode("utf-8").splitlines()
        else:
            raise TypeError("file_content must be a file path or a file-like object.")
        
        self.data = {}
        self._parse_file()
        self._validate_keys()
    
    def _parse_file(self):
        bohr2ang = 0.529177
        lines = self.lines

        for i in range(len(lines)):
            if "siesta: Atomic coordinates (Bohr) and species" in lines[i]:
                pos_cart = []; ion_types = [] ; j = i+1
                while len(lines[j])>1:
                    linesplit = lines[j].split() ; j += 1
                    pos_cart.append(linesplit[1:4]) ; ion_types.append(linesplit[4])
                elements, nions = np.unique(ion_types, return_counts=True)
                self.data["nions"] = np.array(nions) ; natoms = sum(np.array(nions))
                self.data["pos_cart"] = np.array(pos_cart).astype("float")*bohr2ang
            elif "%block ChemicalSpeciesLabel" in lines[i]:
                species = []
                j = i+1
                while "%endblock ChemicalSpeciesLabel" not in lines[j]:
                    linesplit = lines[j].split() ; species.append(linesplit[2])
                    j += 1
                self.data["species"] = species
            elif "outcell: Unit cell vectors (Ang):" in lines[i]:
                lattice = []
                for j in range(i+1, i+4):
                    lattice.append(list(map(float, lines[j].split()))[:3])
                self.data["lattice"] = np.array(lattice)
            elif "Mulliken Atomic Populations:" in lines[i]:
                if len(lines[i+2].split()) == 5:
                    magmom = []
                    for j in range(i+2, i+2+natoms):
                        linesplit = lines[j].split() ; magmom.append([0, 0, linesplit[3]])
                elif len(lines[i+2].split()) == 8:
                    magmom = []
                    for j in range(i+2, i+2+natoms):
                        linesplit = lines[j].split() ; magmom.append(linesplit[4:7])
                else:
                    print("ERROR: Incorret number of spin components !")
                self.data["magmom"] = np.array(magmom).astype("float")

    def _validate_keys(self):
        required_keys = ["lattice", "pos_cart", "nions", "species", "magmom"]
        missing_keys = [key for key in required_keys if key not in self.data]
        if missing_keys:
            raise ValueError(f"Missing required keys in data: {missing_keys}")
    
    def get(self, key, default=None):
        return self.data.get(key, default)

