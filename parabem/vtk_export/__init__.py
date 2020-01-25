import os

import parabem
import numpy

class CaseToVTK():
    def __init__(self, case, _dir, suffix=None):
        """automatic export of a case"""
        self.case = case
        self.make_sure_path_exists(_dir)
        self._dir = _dir
        self.suffix = ""
        if suffix:
            self.suffix = "_" + suffix

    def write_panels(self, cp=True, vel=True, pot=True, data_type="cell"):
        """type=cell/point"""
        with open(self._dir + "/panels" + self.suffix + ".vtk", "w") as _file:
            writer = VtkWriter()
            writer.unstructed_grid(_file, "panels")
            writer.points(_file, self.vertices)
            writer.flat_cells(_file, self.panels)
            if data_type is not "point":
                if cp:
                    writer.data(_file, self.pan_cp, "cp")
                if pot:
                    writer.data(_file, self.pan_pot, "potential")
                if vel:
                    writer.data(_file, self.pan_vel, "velocity", _type="VECTORS")
            else:
                if cp:
                    writer.data(_file, self.vert_cp, "cp", data_type="POINT_DATA")
                if pot:
                    writer.data(_file, self.vert_pot, "potential", data_type="POINT_DATA")
                if vel:
                    writer.data(_file, self.vert_vel, "velocity", _type="VECTORS", data_type="POINT_DATA")

    def write_wake_panels(self, pot=True, normal=True):
        with open(self._dir + "/wakepanels" + self.suffix + ".vtk", "w") as _file:
            writer = VtkWriter()
            writer.unstructed_grid(_file, "wakepanels")
            writer.points(_file, self.wake_vertices)
            writer.flat_cells(_file, self.wake_panels)
            if normal:
                writer.data(_file, [i.n for i in self.case.wake_panels], "normals", _type="VECTORS")
            if pot:
                writer.data(_file, self.wake_pot, "potential")

    def write_field(self, x_seq, y_seq, z_seq, vel=True, pot=True, perturbation=False):
        dim = [x_seq[-1], y_seq[-1], z_seq[-1]]
        x_grid = numpy.linspace(*x_seq)
        y_grid = numpy.linspace(*y_seq)
        z_grid = numpy.linspace(*z_seq)
        vertices = [parabem.PanelVector3(x, y, z) for z in z_grid for y in y_grid for x in x_grid]
        for vert in vertices:
            if pot:
                self.case.off_body_potential(vert)
            if vel or perturbation:
                self.case.off_body_velocity(vert)
        with open(self._dir + "/field.vtk", "w") as _file:
            writer = VtkWriter()
            writer.structed_grid(_file, "field_data", dim)
            writer.points(_file, vertices)
            if pot:
                _potential = [vert.potential for vert in vertices]
                writer.data(_file, _potential, "potential", data_type="POINT_DATA")
            if vel:
                _velocity = [vert.velocity for vert in vertices]
                writer.data(_file, _velocity, "velocity", _type="VECTORS", data_type="POINT_DATA")
            if perturbation:
                _velocity = [vert.velocity.x - self.case.v_inf.x for vert in vertices]
                writer.data(_file, _velocity, "velocity", data_type="POINT_DATA")

    def write_stream_lines(self, start_points=[[0, 0, 0]], interval=0.01, numpoints=100):
        vertices, line_numbers = self._stream_lines(start_points, interval, numpoints)
        with open(self._dir + "/stream_lines.vtk", "w") as _file:
            writer = VtkWriter()
            writer.unstructed_grid(_file, "flow_path")
            writer.points(_file, vertices)
            writer.lines(_file, line_numbers)

    def write_body_stream(self, start_panels, num_panels=20):
        vertices, line_numbers = self._body_stream_lines(start_panels, num_panels)
        with open(self._dir + "/body_stream_lines.vtk", "w") as _file:
            writer = VtkWriter()
            writer.unstructed_grid(_file, "flow_path")
            writer.points(_file, vertices)
            writer.lines(_file, line_numbers)

    @staticmethod
    def make_sure_path_exists(path):
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    @property
    def vertices(self):
        """set the vertex numbers"""
        return self.case.vertices

    @property
    def wake_vertices(self):
        return [vert for pan in self.case.wake_panels for vert in pan.points]

    @property
    def panels(self):
        return [[vert.nr for vert in pan.points] for pan in self.case.panels]

    @property
    def wake_panels(self):
        pans_nr = []
        i = 0
        for pan in self.case.wake_panels:
            verts_nr = []
            for vert in pan.points:
                verts_nr.append(i)
                i += 1
            pans_nr.append(verts_nr)
        return(pans_nr)

    @property
    def pan_vel(self):
        return [pan.velocity for pan in self.case.panels]

    @property
    def vert_vel(self):
        return [vert.velocity for vert in self.case.vertices]

    @property
    def pan_cp(self):
        return [pan.cp for pan in self.case.panels]

    @property
    def vert_cp(self):
        return [vert.cp for vert in self.case.vertices]

    @property
    def wake_pot(self):
        return [pan.potential for pan in self.case.wake_panels]

    @property
    def pan_pot(self):
        return [pan.potential for pan in self.case.panels]

    @property
    def vert_pot(self):
        return [vert.potential for vert in self.case.vertices]

    def _stream_lines(self, start_points=[[0, 0, 0]], interval=0.01, numpoints=100):
        if not isinstance(start_points[0], parabem.Vector3):
            start_points = list(map(parabem.Vector3, start_points))
        print("COMPUTE STREAMLINES")
        flow_paths = [self.case.flow_path(i, interval, numpoints) for i in start_points]
        vertices = []
        fpw = []
        i = 0
        for fp in flow_paths:
            fpwi = []
            for vertex in fp:
                vertices.append(vertex)
                fpwi.append(i)
                i += 1
            fpw.append(fpwi)
        return vertices, fpw

    def _body_stream_lines(self, start_panels, numpanels=100):
        flow_paths = [self.case.body_flow_path(i, numpanels) for i in start_panels]
        vertices = []
        fpw = []
        i = 0
        for fp in flow_paths:
            fpwi = []
            for vertex in fp:
                vertices.append(vertex)
                fpwi.append(i)
                i += 1
            fpw.append(fpwi)
        return vertices, fpw


class VtkWriter():
    def __init__(self):
        self._data_set = False

    def structed_grid(self, _file, name, dim):
        _file.write("# vtk DataFile Version 3.0\n")
        _file.write(name + " \n")
        _file.write("ASCII\n")
        _file.write("DATASET STRUCTURED_GRID\n")
        _file.write("DIMENSIONS ")
        for dim_i in dim:
            _file.write(str(dim_i) + " ")


    def unstructed_grid(self, _file, name):
        _file.write("# vtk DataFile Version 3.0\n")
        _file.write(name + " \n")
        _file.write("ASCII\n")
        _file.write("DATASET UNSTRUCTURED_GRID\n")

    def points(self, _file, vertices):
        _file.write("POINTS " + str(len(vertices)) + " float\n")
        for point in vertices:
            v = parabem.Vector3(point)
            _file.write(str(v[0]) + " " +
                        str(v[1]) + " " +
                        str(v[2]) + "\n")

    def flat_cells(self, _file, cells):
        _file.write("\nCELLS " + str(len(cells)) + " ")
        n = 0
        for i in cells:
            n += len(i) + 1
        _file.write(str(n) + "\n")
        for i in cells:
            _file.write(str(len(i)))
            for j in i:
                _file.write(" " + str(j))
            _file.write("\n")

        _file.write("\nCELL_TYPES " + str(len(cells)) + "\n")
        for j, i in enumerate(cells):
            if j % 5 == 0:
                _file.write("\n")
            if len(i) == 3:
                _file.write("5 ")
            elif len(i) == 4:
                _file.write("9 ")
            else:
                _file.write("7 ")
        _file.write("\n")


    def data(self, _file, data, name="data", _type="SCALARS", data_type="CELL_DATA"):
        if not self._data_set:
            _file.write(data_type +" " + str(len(data)) + "\n")
            self._data_set = True

        if _type is "SCALARS":
            _file.write("SCALARS " + name + " float\n")
            _file.write("LOOKUP_TABLE default")
            for j, value in enumerate(data):
                if j % 5 == 0:
                    _file.write("\n")
                _file.write(str(value) + " ")
            _file.write("\n\n")

        if _type is "VECTORS":
            _file.write("VECTORS " + name + " float\n")
            _file.write("\n")
            for point in data:
                v = parabem.Vector3(point)
                _file.write(
                    str(v[0]) + " " +
                    str(v[1]) + " " +
                    str(v[2]) + "\n")
            _file.write("\n\n")


    def lines(self, _file, lines):
        """lines = [[0,1,2,3,4,5], [6,7,8,9,10]...]"""
        n_lines = sum(map(len, lines)) - len(lines)
        _file.write("\nCELLS " + str(n_lines) + " " + str(n_lines * 3) + "\n")
        for line in lines:
            for j, k in enumerate(line[:-1]):
                _file.write("2 ")
                _file.write(str(k) + " " + str(line[j + 1]))
                _file.write("\n")

        _file.write("\nCELL_TYPES " + str(n_lines) + "\n")
        for k, j in enumerate(range(n_lines)):
            if j % 5 == 0:
                _file.write("\n")
            _file.write("3 ")
        _file.write("\n")