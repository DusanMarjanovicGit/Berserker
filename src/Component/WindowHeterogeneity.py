class WindowHeterogeneity:
    def __init__(self,e_x,e_y,e_z,macroRegionName:str):
        self.e_x = list(sorted(e_x))
        self.e_y = list(sorted(e_y))
        self.e_z = list(sorted(e_z))
        self.macroRegionName = macroRegionName

    def overlaps(self, other):
        """
        Check if this WindowHeterogeneity overlaps with another in 3D space.

        Parameters:
        - other (WindowHeterogeneity): The other instance to compare against.

        Returns:
        - bool: True if the two instances overlap, False otherwise.
        """
        x_overlap = not (self.e_x[1] <= other.e_x[0] or self.e_x[0] >= other.e_x[1])
        y_overlap = not (self.e_y[1] <= other.e_y[0] or self.e_y[0] >= other.e_y[1])
        z_overlap = not (self.e_z[1] <= other.e_z[0] or self.e_z[0] >= other.e_z[1])
        return x_overlap and y_overlap and z_overlap

    def is_included(self, other):
        """Check if this window is fully included in another."""
        return (
            self.e_x[0] >= other.e_x[0] and self.e_x[1] <= other.e_x[1] and
            self.e_y[0] >= other.e_y[0] and self.e_y[1] <= other.e_y[1] and
            self.e_z[0] >= other.e_z[0] and self.e_z[1] <= other.e_z[1]
        )
    def size(self):
        """Calculate the volume of the heterogeneity window."""
        return (self.e_x[1] - self.e_x[0]) * (self.e_y[1] - self.e_y[0]) * (self.e_z[1] - self.e_z[0])
    def __eq__(self, other):
        return (
            isinstance(other, WindowHeterogeneity) and
            self.e_x == other.e_x and
            self.e_y == other.e_y and
            self.e_z == other.e_z and
            self.macroRegionName == other.macroRegionName
        )

    def __hash__(self):
        return hash((tuple(self.e_x), tuple(self.e_y), tuple(self.e_z), self.macroRegionName))

    @staticmethod

    def build_hierarchy(windows):
        """
        Build a hierarchical structure of windows.

        Parameters:
        - windows (list of WindowHeterogeneity): The list of windows to process.

        Returns:
        - list of WindowHeterogeneity: A list of root windows, each with a hierarchy of children.
        """
        # Sort by size (largest first)
        windows = sorted(windows, key=lambda w: w.size(), reverse=True)

        def find_children(parent, candidates):
            """Recursively find children for a given parent window."""
            children = [w for w in candidates if w.is_included(parent) and w != parent]
            for child in children:
                remaining_candidates = [w for w in candidates if w != child]
                child.children = find_children(child, remaining_candidates)
            return children

        roots = []
        while windows:
            root = windows[0]
            windows = [w for w in windows if w != root]
            root.children = find_children(root, windows)
            roots.append(root)

        return roots

    @staticmethod
    def resolve_overlaps(windows):
        """
        Resolve overlaps among a list of WindowHeterogeneity instances by prioritizing larger windows.

        Parameters:
        - windows (list of WindowHeterogeneity): The list of windows to process.

        Returns:
        - list of WindowHeterogeneity: A filtered list of windows with resolved overlaps.
        """
        # Sort windows by size (largest first)
        windows = sorted(windows, key=lambda w: w.size(), reverse=True)

        resolved = []
        for window in windows:
            if all(not window.overlaps(w) for w in resolved):
                resolved.append(window)

        return resolve

    def __repr__(self) : return f"Xbounds: {self.e_x} , YBounds = {self.e_y} , ZBounds = {self.e_z} , xs = {self.macroRegionName}"


