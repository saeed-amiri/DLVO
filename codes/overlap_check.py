
class QuadNode:
    """Represents a node in the Quadtree.
    WiKi:
        A quadtree is a tree data structure in which each internal node
        has exactly four children. Quadtrees are the two-dimensional
        analog of octrees and are most often used to partition a
        two-dimensional space by recursively subdividing it into four
        quadrants or regions. The data associated with a leaf cell
        varies by application, but the leaf cell represents a "unit of
        interesting spatial information".
        The subdivided regions may be square or rectangular, or may
        have arbitrary shapes.
        All forms of quadtrees share some common features:
            - They decompose space into adaptable cells
            - Each cell (or bucket) has a maximum capacity. When
              maximum capacity is reached, the bucket splits
            - The tree directory follows the spatial decomposition of
              the quadtree.
    
    Setting up a complete quadtree structure will involve several steps.
    Quadtrees are a type of tree data structure, typically used in 2D
    spatial partitioning. Let's break down the implementation into
    parts:
        - Quadtree Node: This will represent a node in the quadtree.
        - Quadtree Structure: This will handle operations such as
          insertion, querying, and splitting.
        - Particle Placement: This will leverage the quadtree structure
          to determine if a particle overlaps with others.
    """

    def __init__(self,
                 pos_x,
                 pos_y,
                 width,
                 height,
                 capacity
                 ) -> None:
        """
        Initializes a new quadtree node.

        Args:
            x (float): The x-coordinate of the node's center.
            y (float): The y-coordinate of the node's center.
            width (float): Half of the node's width.
            height (float): Half of the node's height.
            capacity (int): Maximum number of points the node can hold
                            before subdividing.

        Attributes:
            points (list): List of particles within this node.
            orthwest, northeast, southwest, and southeast (QuadNode):
                Child nodes, representing northwest, northeast,
                southwest, and southeast quadrants respectively.
                Initialized as None.
        """
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.width = width
        self.height = height
        self.capacity = capacity
        self.particles = []
        self.divided = False
        self.northwest = None
        self.northeast = None
        self.southwest = None
        self.southeast = None

    def insert(self,
               particle
               ) -> bool:
        """
        Inserts a particle into the quadtree. 

        If the node exceeds its capacity after the insertion, it
        subdivides into four child nodes.
        
        Args:
            particle (Particle): The particle to be inserted.

        Returns:
            bool: True if the point was inserted successfully,
                  False otherwise.
        """
        if not self._contains(particle):
            return False

        if len(self.particles) < self.capacity:
            self.particles.append(particle)
            return True

        if not self.divided:
            self._subdivide()

        if (self.northwest.insert(particle) or
            self.northeast.insert(particle) or
            self.southwest.insert(particle) or
            self.southeast.insert(particle)):
            return True

        return False

    def _contains(self,
                  particle
                  ) -> bool:
        """
        Checks if the given particle lies within the boundaries of
        this node.

        Args:
            particle (Particle): The particle to check.

        Returns:
            bool: True if the particle is inside the node, False
            otherwise.
        """
        return (
            self.pos_x - self.width < particle.position[0] <= 
            self.pos_x + self.width) and (
            self.pos_y - self.height < particle.position[1] <=
            self.pos_y + self.height)

    def _subdivide(self) -> None:
        northwest = QuadNode(self.pos_x - self.width / 2,
                             self.pos_y - self.height / 2,
                             self.width / 2,
                             self.height / 2,
                             self.capacity)
        northeast = QuadNode(self.pos_x + self.width / 2,
                             self.pos_y - self.height / 2,
                             self.width / 2,
                             self.height / 2,
                             self.capacity)
        southwest = QuadNode(self.pos_x - self.width / 2,
                             self.pos_y + self.height / 2,
                             self.width / 2,
                             self.height / 2,
                             self.capacity)
        southeast = QuadNode(self.pos_x + self.width / 2,
                             self.pos_y + self.height / 2,
                             self.width / 2,
                             self.height / 2,
                             self.capacity)

        self.northwest = northwest
        self.northeast = northeast
        self.southwest = southwest
        self.southeast = southeast

        self.divided = True
