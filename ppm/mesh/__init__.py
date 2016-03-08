from __future__ import division
import ppm

# exdent with other mesh formats (eg. stl)

class mesh_object(object):
    def __init__(self):
        self.vertices = []
        self.panels = []
        self.trailing_edges = []        

    @classmethod
    def from_OBJ(cls, path):

        def load_obj_file(path):
            trailing_edges = []
            panels = []
            vertices = []
            with open(path, "r") as inp:
                text_list = inp.readlines()
            for line in text_list:
                chars = line.split()
                if chars[0] == "v":
                    vertices.append([float(number) for number in chars[1:]])
                elif chars[0] == "f":
                    panels.append([int(number) - 1 for number in chars[1:]])
                elif chars[0] == "l":
                    trailing_edges.append(
                        [int(number) - 1 for number in chars[1:]])
            return vertices, panels, trailing_edges

        def sort_and_order(panels, vertices):
            # deleting all vertices not used by the mesh
            # overwrites the panels and return the new vertices list
            orderd_vert_nr = list(set([vert for pol in panels for vert in pol]))
            vertices_set = [vertices[i] for i in orderd_vert_nr]
            mesh_dict = dict(zip(orderd_vert_nr, range(len(orderd_vert_nr))))
            panels = [[mesh_dict[vert] for vert in pol] for pol in panels]
            return panels, vertices_set

        def sort_wake(trailing_edges):
            # makes a orderd list of the trailing_edge
            w_temp = trailing_edges[0]
            trailing_edges.pop(0)
            i = 0
            j = 0  # stop if there is an error in the data
            while len(trailing_edges) and j < 10000:
                if trailing_edges[i][0] == w_temp[0]:
                    w_temp = [trailing_edges[i][1]] + w_temp
                    trailing_edges.pop(i)
                elif trailing_edges[i][1] == w_temp[0]:
                    w_temp = [trailing_edges[i][0]] + w_temp
                    trailing_edges.pop(i)
                elif trailing_edges[i][0] == w_temp[-1]:
                    w_temp = w_temp + [trailing_edges[i][1]]
                    trailing_edges.pop(i)
                elif trailing_edges[i][1] == w_temp[-1]:
                    w_temp = w_temp + [trailing_edges[i][0]]
                    trailing_edges.pop(i)
                i = int((i + 1) * (i < (len(trailing_edges) - 1)))
                j += 1
            # abbort if j == 10000: this means that the trailing_edge is not a single line
            # this has to be exdented to wakes consisting of more than one wake edge !!!
            return w_temp

        def wake_from_mesh(trailing_edges, vertices, vertices_set):
            # replacing vertices with mesh vertices
            norm = lambda x, y: ((x[0] - y[0])**2 +
                                 (x[1] - y[1])**2 +
                                 (x[2] - y[2])**2)**(0.5)
            vert_dict = dict(zip(range(len(vertices)), vertices))
            max_size = 0.0001
            wake_dict = {}
            for i in trailing_edges:
                v1 = vert_dict[i]
                for j, v2 in enumerate(vertices_set):
                    if norm(v1, v2) < max_size:
                        wake_dict[i] = j
            return wake_dict

        mesh = cls()
        vertices, panels, trailing_edges = load_obj_file(path)
        panels, vertices_set = sort_and_order(panels, vertices)
        mesh.vertices = [ppm.PanelVector3(*vertex) for vertex in vertices_set]
        mesh.panels = [ppm.Panel3([mesh.vertices[nr] for nr in pol]) for pol in panels]
        if len(trailing_edges) > 0:
            trailing_edges = sort_wake(trailing_edges)
            wake_dict = wake_from_mesh(trailing_edges, vertices, vertices_set)
            mesh.trailing_edges = [mesh.vertices[wake_dict[i]] for i in trailing_edges]
        return mesh
