import sys
import io
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
global matrix_mod_rot_1, ma_f_1_rot_1, ma_f_2_rot_1, header_change_0, header_change_1, header_change_2, temp_input_data, coef_val_2, max_value, f_1, f_2, resp_f, \
    f_1max, f_2max, my_results, results_buffer_1, results_buffer, Results, r_data_input, Results_2, Results, select, n_replic, n_rows, designs_list, ma_f_1, matrix_f_1, matrix_f_2, ma_f_1, rand_index, matrix_mod, matrix_f_1_1, matrix_f_2_1, ma_f_1_1, rand_index_1, matrix_mod_1, matrix_f_2_m, f_1_max, f_2_max
# Very important, Define all the globals! :))
# Name of factors and response
f_1_name = ""
f_2_name = ""
resp_name = ""
f_1_name_temp = "FACTOR 1"
f_2_name_temp = "FACTOR 2"
resp_name_temp = "RESPONSE"
# --------------------------------------------------------------------------------------------------------------------
# External Functions

# Extracting data from buffer
def extracting_data_from_data_frame():
    global f_1, f_2, resp_f, temp_input_data, r_data_input
    # Extracting column data from Data Frame stored in buffer (temp_input_data)
    r_data_input = pd.read_csv(io.StringIO(temp_input_data.getvalue()))
    f_1 = r_data_input["0"].explode()
    f_2 = r_data_input["1"].explode()
    resp_f = r_data_input["2"].explode()
    return
# Processing the input data
def processing():
    global temp_input_data, f_1, f_2, resp_f
    # Conditional Analysis, the plot will be done only when number of experimental points is higher than 3
    if len(test) > 3:
        # Extracting column data from Data Frame stored in buffer (temp_input_data)
        extracting_data_from_data_frame()
    return
def fact_names():
    global f_1_name, f_2_name, resp_name, f_1_name_temp, f_2_name_temp, resp_name_temp
    if f_1_name_temp == "":
        f_1_name_temp = "FACTOR 1"
    if f_1_name == "":
        f_1_name = f_1_name_temp
    if f_2_name_temp == "":
        f_2_name_temp = "FACTOR 2"
    if f_2_name == "":
        f_2_name = f_2_name_temp
    if resp_name_temp == "":
        resp_name_temp = "RESPONSE"
    if resp_name == "":
        resp_name = resp_name_temp
    return

# Classes
class initial_window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.gd = None
        self.setGeometry(20, 20, 800, 600)
        self.setFixedSize(720, 550)
        self.setWindowTitle("M-Design: A project for Design of Experiments, ver 1.0, by Armando Hernandez (2022)")
        self.insert_table = input_table(100, 3)
        self.setCentralWidget(self.insert_table)
        # Adding a Menu Bar with submenus
        my_bar = QToolBar()
        my_bar.setIconSize(QSize(16, 16))
        font_bar = QFont()
        font_bar.setPointSize(15)
        open_file = QAction(QIcon("/PATH/blue-folder-horizontal-open.png"), "Open File", self)  # Open file
        open_file.setFont(font_bar)
        open_file.setChecked(True)
        open_file.triggered.connect(self.insert_table.open_file)
        save_file = QAction(QIcon("/PATH/disk-black.png"), "Save File as CSV", self)  # Save File
        save_file.setFont(font_bar)
        save_file.setChecked(True)
        save_file.triggered.connect(self.insert_table.save_data)
        define_factors = QAction(QIcon("/PATH/flask.png"), "Define Factors", self)  # Define Factors
        define_factors.setFont(font_bar)
        define_factors.setChecked(True)
        define_factors.triggered.connect(self.insert_table.factors_names_connection)
        gen_design = QAction(QIcon("/PATH/generate design-2.png"), "Generate Design", self)  # Generate Design
        gen_design.setFont(font_bar)
        gen_design.setChecked(True)
        gen_design.triggered.connect(self.gene_design)
        scatter_plot = QAction(QIcon("/PATH/3d-modeling.png"), "Scatter Plot", self)  # Scatter Plot
        scatter_plot.setFont(font_bar)
        scatter_plot.setCheckable(True)
        scatter_plot.triggered.connect(self.insert_table.data_in_table_for_plot_scattered)
        linear_fit = QAction(QIcon("/PATH/linear fitting.png"), "Fit to Linear", self)  # Linear Fit
        linear_fit.setFont(font_bar)
        linear_fit.setChecked(True)
        linear_fit.triggered.connect(self.insert_table.data_in_table_for_fit_first_order)
        quadratic_fit = QAction(QIcon("/PATH/quadratic fitting.png"), "Fit to quadratic", self)  # Quadratic Fit
        quadratic_fit.setFont(font_bar)
        quadratic_fit.setChecked(True)
        quadratic_fit.triggered.connect(self.insert_table.data_in_table_for_fit_quadratic)
        plot_surface_1 = QAction(QIcon("/PATH/Plot-surface1"), "Plot Surface-1", self)  # Plot Surface-1
        plot_surface_1.setFont(font_bar)
        plot_surface_1.setChecked(True)
        plot_surface_1.triggered.connect(self.insert_table.data_in_table_for_surface_first_order)
        plot_surface_2 = QAction(QIcon("/PATH/Surface2.png"), "Plot Surface-2", self)  # Plot Surface-2
        plot_surface_2.setFont(font_bar)
        plot_surface_2.setChecked(True)
        plot_surface_2.triggered.connect(self.insert_table.data_in_table_for_surface_quadratic)
        clear_cells = QAction("Clear Cells", self)  # Clear Cells (still the connect button is missing)
        clear_cells.setFont(font_bar)
        clear_cells.setChecked(True)
        clear_cells.triggered.connect(self.insert_table.clear_cells_input_table)
        clear_all = QAction("Clear All", self)  # Clear All (still the connect button is missing)
        clear_all.setFont(font_bar)
        clear_all.setChecked(True)
        clear_all.triggered.connect(self.insert_table.clear_all)
        quit_from = QAction("Quit", self)  # Quit (still the connect button is missing)
        quit_from.setFont(font_bar)
        quit_from.setChecked(True)
        quit_from.triggered.connect(self.quit_0)
        icon_bar = [open_file, save_file, define_factors, gen_design, scatter_plot, linear_fit, quadratic_fit, plot_surface_1, plot_surface_2, clear_cells, clear_all, quit_from]
        self.addToolBar(my_bar)
        my_bar.addActions(icon_bar)
        # Help Option
        #button_2 = QPushButton(self)
        #button_2.setGeometry(600, 230, 80, 30)
        #button_2.setFont(font_0)
        #button_2.setText("Help")
        #button_2.clicked.connect(self.help_description)
        self.show()
# Def Quit
    def quit_0(self):
        self.insert_table.clear_all()
        self.insert_table.close()
        self.close()
# Def gene_design
    def gene_design(self):
        if self.gd is None:
            self.gd = generate_design()
            self.gd.show()
        elif generate_design():
            self.gd = generate_design()
            self.gd.show()
        else:
            self.gd.close()
            self.gd = None
# Not Enough Data!
class not_enough_data(QDialog):
    def __init__(self):
        super().__init__()
        self.setGeometry(500, 300, 450, 120)
        self.setFixedSize(450, 120)
        font = QFont()
        font.setPointSize(15)
        font.setBold(True)
        self.setFont(font)
        n_e_data = QLabel(self)
        n_e_data.setGeometry(40, 20, 400, 30)
        n_e_data.setText("Not enough DATA to process,  please add more DATA!")
        ok_button = QPushButton(self)
        ok_button.setGeometry(170, 60, 60, 30)
        ok_button.setText("OK")
        ok_button.clicked.connect(self.close)
        self.show()
# Generate Design
class generate_design(QWidget):
    def __init__(self):
        global select, n_replic, n_rows, designs_list
        n_replic = 2
        super().__init__()
        style = QVBoxLayout()
        self.design = None
        self.setWindowTitle("Design Generation")
        self.setGeometry(800, 20, 300, 280)
        self.setFixedSize(300, 280)
        font = QFont()
        font.setPointSize(14)
        font.setBold(True)
        self.type_design = QLabel(self)
        self.replicates = QLabel(self)
        self.type_design.setFont(font)
        self.type_design.setGeometry(80, 50, 300, 20)
        self.type_design.setText("Choose your design")
        self.replicates.setFont(font)
        self.replicates.setGeometry(80, 140, 300, 20)
        self.replicates.setText("Number of replicates")
        self.box = QComboBox()
        self.box.setFont(font)
        designs_list = ["Full Factorial (Two levels)", "Full Factorial (Three levels)", "Central Composite Design"]
        self.box.addItems(designs_list)
        self.n_replicates = QComboBox()
        self.n_replicates.setFont(font)
        self.n_replicates.addItems(["2", "3", "4", "5"])
        self.n_replicates.currentIndexChanged.connect(self.replic)
        self.cont = QPushButton(self)
        self.cont.setGeometry(100, 230, 100, 30)
        self.cont.setText("OK")
        self.cont.clicked.connect(self.design_selection)
        style.addWidget(self.box)
        style.addWidget(self.n_replicates)
        self.setLayout(style)
        self.show()

    # Number of Replicates Selection
    def replic(self):
        global n_replic
        n_replic = int(self.n_replicates.currentText())
        return

    # Design selection
    def design_selection(self):
        global select, designs_list, n_rows
        select = self.box.currentText()
        for x in range(0, len(designs_list)):
            if select == designs_list[x]:
                self.design = design_matrix()
                self.design.show()
        return
class table_design_matrix(QTableWidget):
    def __init__(self):
        global matrix_mod_rot_1, ma_f_1_rot_1, ma_f_2_rot_1, designs_list, n_replic, n_rows, f_1_name, f_2_name, resp_name, header_change_0, header_change_1, header_change_2, select, matrix_f_1, matrix_f_2, ma_f_1, rand_index, matrix_mod, matrix_f_1_1, matrix_f_2_1, ma_f_1_1, rand_index_1, matrix_mod_1, matrix_f_2_m_1, matrix_f_2_m
        fact_names()
        self.analysis()
        super().__init__(n_rows, 3)
        self.iw = None
        self.bm = None
        self.fn = None
        self.w = None
        self.w1 = None
        self.ps = None
        self.fo = None
        self.fop = None
        self.qo = None
        self.qop = None
        self.ne = None
        self.setWindowTitle("Design Matrix")
        self.setGeometry(450, 200, 450, 300)
        self.setFixedSize(400, 350)
        self.setColumnWidth(0, 100)
        self.setColumnWidth(1, 100)
        self.setColumnWidth(2, 100)
        self.init_ui()
        # Change of Name Factor 1
        f_1_name, header_change_0 = QInputDialog.getText(self, "Factor 1 Name", 'FACTOR 1:', QLineEdit.Normal, f_1_name)
        if f_1_name == "":
            f_1_name = f_1_name_temp
        if header_change_0:
            self.setHorizontalHeaderItem(0, QTableWidgetItem(f_1_name))
        # Change of Name Factor 2
        f_2_name_temp = f_2_name
        f_2_name, header_change_1 = QInputDialog.getText(self, "Factor 2 Name", 'FACTOR 2:', QLineEdit.Normal, f_2_name)
        if f_2_name == "":
            f_2_name = f_2_name_temp
        if header_change_1:
            self.setHorizontalHeaderItem(1, QTableWidgetItem(f_2_name))
        # Change of Name Response
        resp_name, header_change_2 = QInputDialog.getText(self, "Response Name", 'RESPONSE:', QLineEdit.Normal, resp_name)
        if resp_name == "":
            resp_name = resp_name_temp
        if header_change_2:
            self.setHorizontalHeaderItem(2, QTableWidgetItem(resp_name))
        # Inserting items in table (combinations of different levels)
        # Generating the levels combination at random
        if select == designs_list[0]:
            self.two_two_calculations()
            for d in range(0, n_rows):
                ma_f_1_1.append(QTableWidgetItem(str(matrix_mod_1[d][0])))
                self.setItem(d, 0, ma_f_1_1[d])
                ma_f_1_1[d].setTextAlignment(Qt.AlignCenter)
                matrix_f_2_m_1.append(QTableWidgetItem(str(matrix_mod_1[d][1])))
                self.setItem(d, 1, matrix_f_2_m_1[d])
                matrix_f_2_m_1[d].setTextAlignment(Qt.AlignCenter)
                ma_f_1_1[d].setFlags(ma_f_1_1[d].flags() ^ Qt.ItemIsEditable)
                matrix_f_2_m_1[d].setFlags(matrix_f_2_m_1[d].flags() ^ Qt.ItemIsEditable)
                self.setVerticalHeaderItem(d, QTableWidgetItem("experiment " + str(d + 1)))
            return
        if select == designs_list[1]:
            self.three_two_calculations()
            for d in range(0, n_rows):
                ma_f_1.append(QTableWidgetItem(str(matrix_mod[d][0])))
                self.setItem(d, 0, ma_f_1[d])
                ma_f_1[d].setTextAlignment(Qt.AlignCenter)
                matrix_f_2_m.append(QTableWidgetItem(str(matrix_mod[d][1])))
                self.setItem(d, 1, matrix_f_2_m[d])
                matrix_f_2_m[d].setTextAlignment(Qt.AlignCenter)
                ma_f_1[d].setFlags(ma_f_1[d].flags() ^ Qt.ItemIsEditable)
                matrix_f_2_m[d].setFlags(matrix_f_2_m[d].flags() ^ Qt.ItemIsEditable)
                self.setVerticalHeaderItem(d, QTableWidgetItem("experiment " + str(d + 1)))
        if select == designs_list[2]:
            self.rotational_1_calculations()
            for d in range(0, n_rows):
                ma_f_1_rot_1.append(QTableWidgetItem(str(matrix_mod_rot_1[d][0])))
                self.setItem(d, 0, ma_f_1_rot_1[d])
                ma_f_1_rot_1[d].setTextAlignment(Qt.AlignCenter)
                ma_f_2_rot_1.append(QTableWidgetItem(str(matrix_mod_rot_1[d][1])))
                self.setItem(d, 1, ma_f_2_rot_1[d])
                ma_f_2_rot_1[d].setTextAlignment(Qt.AlignCenter)
                ma_f_1_rot_1[d].setFlags(ma_f_1_rot_1[d].flags() ^ Qt.ItemIsEditable)
                ma_f_2_rot_1[d].setFlags(ma_f_2_rot_1[d].flags() ^ Qt.ItemIsEditable)
                self.setVerticalHeaderItem(d, QTableWidgetItem("experiment " + str(d + 1)))
        return

    def analysis(self):
        global n_rows, designs_list, select
        if select == designs_list[0]:
            n_rows = n_replic * 4
            self.two_two_calculations
            return
        if select == designs_list[1]:
            n_rows = n_replic * 9
            self.three_two_calculations
        if select == designs_list[2]:
            n_rows = n_replic * 9
            self.rotational_1_calculations
            return
    # rotational_1
    def rotational_1_calculations(self):
        global matrix_rot_1_1, matrix_rot_1_2,  ma_f_1_rot_1, ma_f_2_rot_1, matrix_mod_rot_1
        matrix_rot_1_1 = [-1, -1.4142, -1.4142, 0, 0, 0, 1.4142, 1.4142, 1]
        matrix_rot_1_2 = [0, -1.4142, 1.4142, 1, -1, 0, -1.4142, 1.4142, 0]
        matrix_rot_1_1_inc = []
        matrix_rot_1_2_inc = []
        ma_f_1_rot_1 = []
        ma_f_2_rot_1 = []
        rand_index = []
        matrix_mod_rot_1 = []
        for i in range(0, n_replic):
            matrix_rot_1_1_inc = matrix_rot_1_1_inc + matrix_rot_1_1
            matrix_rot_1_2_inc = matrix_rot_1_2_inc + matrix_rot_1_2
        for k in range(0, len(matrix_rot_1_2_inc)):
            rand_index = np.random.permutation(k + 1)
        mexp = np.array([matrix_rot_1_1_inc, matrix_rot_1_2_inc])
        mexp = np.transpose(mexp)
        for s in rand_index:
            matrix_mod_rot_1.append(mexp[s])
        return
    # 2^2 Calculations
    def two_two_calculations(self):
        global matrix_f_1_1, matrix_f_2_1, ma_f_1_1, rand_index_1, matrix_mod_1, matrix_f_2_m_1
        # Generating de matrix of both factors
        matrix_f_1_1 = [1, -1, 1, -1]
        matrix_f_2_1 = [1, -1, -1, 1]
        matrix_f_1_1_inc = []
        matrix_f_2_1_inc = []
        ma_f_1_1 = []
        matrix_f_2_m_1 = []
        rand_index_1 = []
        matrix_mod_1 = []
        for i in range(0, n_replic):
            matrix_f_1_1_inc = matrix_f_1_1_inc + matrix_f_1_1
            matrix_f_2_1_inc = matrix_f_2_1_inc + matrix_f_2_1
        # Generating the levels combination at random
        for d in range(0, n_rows):
            rand_index_1 = np.random.permutation(d+1)
        mexp_1 = np.array([matrix_f_1_1_inc, matrix_f_2_1_inc])
        mexp_1 = np.transpose(mexp_1)
        matrix_mod_1 = []
        for s in rand_index_1:
            matrix_mod_1.append(mexp_1[s])
        return
    # 3^2 Calculations
    def three_two_calculations(self):
        global matrix_f_1, matrix_f_2, ma_f_1, rand_index, matrix_mod, n_rows, matrix_f_2_m
        # Generating de matrix of both factors
        matrix_f_1 = []
        matrix_f_2 = []
        ma_f_1 = []
        rand_index = []
        matrix_f_2_m = []
        matrix_mod = []
        for i in range(0, n_rows):
            matrix_f_1.append((-1) ** (i))
            k = 1
            l = 2
            m = 0
            for j in range(0, len(matrix_f_1)):
                if j == k:
                    matrix_f_1[j] = matrix_f_1[k] = 0
                    k = k + 3
                if j == l:
                    matrix_f_1[j] = matrix_f_1[l] = 1
                    l = l + 3
                if j == m:
                    matrix_f_1[j] = matrix_f_1[m] = -1
                    m = m + 3
        # Generating the second paired matrix - combination of levels
        k = 0
        l = 1
        m = 6
        matrix_test = []
        for i in range(0, 9):
            matrix_test.append((-1 ** i))
            if i == k:
                matrix_test[i] = matrix_test[k] = 0
                k = k + 2
            if i == l:
                matrix_test[i] = matrix_test[l] = -1
                l = l + 2
            if i >= m:
                matrix_test[i] = matrix_test[m] = 1
        for j in range(0, n_replic):
            matrix_f_2 = matrix_f_2 + matrix_test
        # Generating the levels combination at random
        for d in range(0, len(matrix_f_2)):
            rand_index = np.random.permutation(d + 1)
        mexp = np.array([matrix_f_1, matrix_f_2])
        mexp = np.transpose(mexp)
        matrix_mod = []
        for s in rand_index:
            matrix_mod.append(mexp[s])
        return
    # Def Init_UI
    def init_ui(self):
        self.cellChanged.connect(self.c_current)
        self.cellEntered.connect(self.c_current)
        self.cellClicked.connect(self.c_current)
        self.cellPressed.connect(self.c_current)
        self.cellDoubleClicked.connect(self.c_current)
        # Context Menu Options
        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        select_option = QAction("Select", self)
        copy_option = QAction("Copy", self)
        paste_option = QAction("Paste", self)
        self.addActions([select_option, copy_option, paste_option])
    # Def Current Cell
    def c_current(self):
        global temp_input_data, result, fact_1, test, n_rows
        result = []
        for col in range(0, 3):
            rows = []
            for row in range(0, n_rows):
                value = self.item(row, col)
                if value:
                    rows.append(value.text())
                    value.setTextAlignment(Qt.AlignCenter)
                    if value.text() == "":
                        rows.clear()
            result.append(rows)
        fact_1 = result[0]
        fact_2 = result[1]
        resp = result[2]
        test = np.array([fact_1, fact_2, resp], dtype=object)
        test = np.transpose(test)
        test_data = pd.DataFrame(test)
        # Saving to buffer
        temp_input_data = io.StringIO()
        test_data.to_csv(temp_input_data)

    # For now insert the functions to calculate. Otherwise, create calculation modules to mport for using in every stage
    #  Def calculations
    def data_in_table_for_plot_scattered(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_scattered_connection()
            return
    # Def Data in table for linear Model with interactions
    def data_in_table_for_fit_first_order(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.first_order_fitting_c()
            return

    # Def Data in table for quadratic Model with interactions
    def data_in_table_for_fit_quadratic(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.second_order_fitting_c()
            return

    # Def Data in table for surface linear Model
    def data_in_table_for_surface_first_order(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_first_order_connection()
            return

    # Def Data in table for surface quadratic Model
    # Plot Surface Decision if no- data entry, then pop-up something...e.g. "You haven't entered data yet..."
    def data_in_table_for_surface_quadratic(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_quadratic_connection()
            return

    # Window to connect the Name of Factors
    def factors_names_connection(self):
        global f_1_name, f_2_name, resp_name, header_change_0, header_change_1, header_change_2, f_1_name_temp, \
            f_2_name_temp, resp_name_temp
        # Change of Name Factor 1
        f_1_name_temp = f_1_name
        f_1_name, header_change_0 = QInputDialog.getText(self, "Factor 1 Name", 'FACTOR 1:', QLineEdit.Normal, f_1_name)
        if f_1_name == "":
            f_1_name = f_1_name_temp
        if header_change_0:
            self.setHorizontalHeaderItem(0, QTableWidgetItem(f_1_name))
        # Change of Name Factor 2
        f_2_name_temp = f_2_name
        f_2_name, header_change_1 = QInputDialog.getText(self, "Factor 2 Name", 'FACTOR 2:', QLineEdit.Normal, f_2_name)
        if f_2_name == "":
            f_2_name = f_2_name_temp
        if header_change_1:
            self.setHorizontalHeaderItem(1, QTableWidgetItem(f_2_name))
        # Change of Name Factor 3
        resp_name_temp = resp_name
        resp_name, header_change_2 = QInputDialog.getText(self, "Response Name", "RESPONSE:", QLineEdit.Normal,
                                                          resp_name)
        if resp_name == "":
            resp_name = resp_name_temp
        if header_change_2:
            self.setHorizontalHeaderItem(2, QTableWidgetItem(resp_name))
        self.c_current()
        return

        # Window to pop up "Not enough data"
    def not_enough_data_connection(self):
        global temp_input_data
        if self.ne is None:
            self.ne = not_enough_data()
            self.ne.show()
        elif not_enough_data():
            self.ne = not_enough_data()
            self.ne.show()
        else:
            self.ne.close()
            self.ne = None

    # Window to show the Scatter Plot
    def plot_scattered_connection(self):
        global temp_input_data
        if self.ps is None:
            self.ps = plot_scattered()
            self.ps.show()
        elif plot_scattered():
            self.ps = plot_scattered()
            self.ps.show()
        else:
            self.ps.close()
            self.ps = None
    # Window to show the First -Order Plot
    def plot_first_order_connection(self):
        if self.fo is None and self.fop is None:
            self.fo = plot_first_order()
            self.fop = plot_first_order_projection()
            self.fo.show()
            self.fop.show()
        elif plot_first_order():
            self.fo = plot_first_order()
            self.fop = plot_first_order_projection()
            self.fo.show()
            self.fop.show()
        else:
            self.fo.close()
            self.fo = None
            self.fop.close()
            self.fop = None

    # Window to show the Second-Order Plot
    def plot_quadratic_connection(self):
        if self.qo is None and self.qop is None:
            self.qo = plot_quadratic()
            self.qop = plot_quadratic_projection()
            self.qo.show()
            self.qop.show()
        elif plot_first_order():
            self.qo = plot_quadratic()
            self.qop = plot_quadratic_projection()
            self.qo.show()
            self.qop.show()
        else:
            self.qo.close()
            self.qop = None
            self.qo.close()
            self.qop = None
    # Window to show the regression results regarding linear fitting
    def first_order_fitting_c(self):
        if self.w is None:
            self.w = first_order_fitting()
            self.w.show()
        elif first_order_fitting():
            self.w = first_order_fitting()
            self.w.show()
        else:
            self.w.close()
            self.w = None

    # Window to show the regression results regarding quadratic fitting
    def second_order_fitting_c(self):
        if self.w1 is None:
            self.w1 = second_order_fitting()
            self.w1.show()
        elif second_order_fitting():
            self.w1 = second_order_fitting()
            self.w1.show()
        else:
            self.w1.close()
            self.w1 = None
# Saving Data as CSV
    def save_data(self):
        global temp_input_data, result, f_1_name, f_2_name, resp_name
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            # Temporary reading data from buffer to save it later
            temp_read = pd.read_csv(io.StringIO(temp_input_data.getvalue()))
            file_path = QFileDialog.getSaveFileName(self, "Save File ", "*.csv")
            file_name = str(file_path)
            str_list = []
            c = 1
            file_name_1 = ""
            for i in range(1, len(file_name)):
                str_list.append(file_name[i])
            for j in range(0, len(str_list)):
                while str_list[c] != "'":
                    file_name_1 = file_name_1 + str_list[c]
                    c = c + 1
            if file_name_1 != "":
                modify_temp_read = temp_read.rename(
                    columns={temp_read.columns[1]: f_1_name, temp_read.columns[2]: f_2_name,
                             temp_read.columns[3]: resp_name})  # Replacing the headers names
                modify_temp_read.to_csv(file_name_1, index=None, columns=[f_1_name, f_2_name, resp_name])
            return


# Two-Two Main Window
class design_matrix(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(450, 200, 400, 300)
        self.setFixedSize(400, 380)
        # Herein try to create different conditions for using the same main window for different alternatives to design.
        # Need a global variable to establish the conditions
        self.setWindowTitle("Design Matrix")
        self.insert_table = table_design_matrix()
        self.setCentralWidget(self.insert_table)
        # Adding a Tool Bar with submenus
        my_bar = QToolBar()
        my_bar.setIconSize(QSize(16, 16))
        font_bar = QFont()
        font_bar.setPointSize(15)
        save_file = QAction(QIcon("/PATH/disk-black.png"), "Save File as CSV", self)  # Save File
        save_file.setFont(font_bar)
        save_file.setChecked(True)
        save_file.triggered.connect(self.insert_table.save_data)
        define_factors = QAction(QIcon("/PATH/flask.png"), "Define Factors", self)  # Define Factors
        define_factors.setFont(font_bar)
        define_factors.setChecked(True)
        define_factors.triggered.connect(self.insert_table.factors_names_connection)
        scatter_plot = QAction(QIcon("/PATH/3d-modeling.png"), "Scatter Plot", self)  # Scatter Plot
        scatter_plot.setFont(font_bar)
        scatter_plot.setCheckable(True)
        scatter_plot.triggered.connect(self.insert_table.data_in_table_for_plot_scattered)
        linear_fit = QAction(QIcon("/PATH/linear fitting.png"), "Fit to Linear", self)  # Linear Fit
        linear_fit.setFont(font_bar)
        linear_fit.setChecked(True)
        linear_fit.triggered.connect(self.insert_table.data_in_table_for_fit_first_order)
        quadratic_fit = QAction(QIcon("/PATH/quadratic fitting.png"), "Fit to quadratic", self)  # Quadratic Fit
        quadratic_fit.setFont(font_bar)
        quadratic_fit.setChecked(True)
        quadratic_fit.triggered.connect(self.insert_table.data_in_table_for_fit_quadratic)
        plot_surface_1 = QAction(QIcon("/PATH/Plot-surface1"), "Plot Surface-1", self)  # Plot Surface-1
        plot_surface_1.setFont(font_bar)
        plot_surface_1.setChecked(True)
        plot_surface_1.triggered.connect(self.insert_table.data_in_table_for_surface_first_order)
        plot_surface_2 = QAction(QIcon("/PATH/Surface2.png"), "Plot Surface-2", self)  # Plot Surface-2
        plot_surface_2.setFont(font_bar)
        plot_surface_2.setChecked(True)
        plot_surface_2.triggered.connect(self.insert_table.data_in_table_for_surface_quadratic)
        quit_from = QAction("Quit", self)
        quit_from.setFont(font_bar)
        quit_from.setChecked(True)
        quit_from.triggered.connect(self.close)
        icon_bar = [save_file, define_factors, scatter_plot, linear_fit, quadratic_fit, plot_surface_1, plot_surface_2, quit_from]
        self.addToolBar(my_bar)
        my_bar.addActions(icon_bar)
        self.show()
# Scattered Plot Canvas
class plot_scattered_canvas(FigureCanvas):
    global temp_input_data, f_1, f_2, resp_f, f_1_name, f_2_name, resp_name, \
        header_change_0, header_change_1, header_change_2
    def __init__(self, parent):
        global temp_input_data, f_1, f_2, resp_f, f_1_name, f_2_name, resp_name, header_change_0, header_change_1, \
            header_change_2
        processing()
        fact_names()
        plt.close(1)
        self.sc_plot = plt.figure(1, figsize=(5, 4))
        super().__init__(self.sc_plot)
        self.setParent(parent)
        self.eje1 = plt.axes(projection="3d")
        self. eje1.scatter(f_1, f_2, resp_f, s=50, c="red")
        self.eje1.set_title("Scatter Plot of Experimental Points")
        self.eje1.set_xlabel(f_1_name)
        self.eje1.set_ylabel(f_2_name)
        self.eje1.set_zlabel(resp_name)
        plt.draw()
        return
# Scattered Plot Main Window
class plot_scattered(QMainWindow):
    def __init__(self):
        super().__init__()
        self.close()
        self.setGeometry(20, 50, 500, 400)
        self.setWindowTitle("Scatter Plot of Experimental Points")
        self.setFixedSize(500, 400)
        plot_scattered_canvas(self)
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setFixedSize(500, 20)
        self.addToolBar(bar)
        save_fig = QAction(QIcon("/PATH/disk-black.png"), "Save Figure", self)
        save_fig.setCheckable(True)
        save_fig.triggered.connect(self.save_scatter_plot)
        bar.addAction(save_fig)
        self.show()
    # Defining Save Figure
    def save_scatter_plot(self):
        scatter_path = QFileDialog.getSaveFileName(self, "Save Scatter Plot ", "*.jpg")
        file_name = str(scatter_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            plt.figure(1).savefig(file_name_1)
        return

# First order Plot Canvas
class plot_first_order_canvas(FigureCanvas):
    # Defining the buffer and variables in the buffer as globals
    global temp_input_data, f_1, f_2, resp_f, const, results_buffer_1
    def __init__(self, parent):
        processing()
        first_order_fitting()
        fact_names()
        plt.close(2)
        my_results = pd.read_csv(io.StringIO(results_buffer_1.getvalue()))
        const = my_results[" Coefficient  "].explode()
        self.sc_plot = plt.figure(2, figsize=(5, 4))
        super().__init__(self.sc_plot)
        self.setParent(parent)
        self.eje1 = plt.axes(projection="3d")
        self. eje1.scatter(f_1, f_2, resp_f, s=50, c="red")
        self.eje1.set_title("Response Surface, Linear Model with Interactions")
        self.eje1.set_xlabel(f_1_name)
        self.eje1.set_ylabel(f_2_name)
        self.eje1.set_zlabel(resp_name)
        s1 = (max(f_1) - min(f_1)) / 100
        s2 = (max(f_2) - min(f_2)) / 100
        x = np.arange(min(f_1), max(f_1), s1)
        y = np.arange(min(f_2), max(f_2), s2)
        X, Y = np.meshgrid(x, y)
        z = const[0] + const[1] * x + const[2] * y + const[3] * x * y
        Z = const[0] + const[1] * X + const[2] * Y + const[3] * X * Y
        self.eje1.plot_surface(X, Y, Z, cmap="winter")
        plt.draw()
        super().__init__(self.sc_plot)
        self.setParent(parent)
# First order Plot Main
class plot_first_order(QMainWindow):
    def __init__(self):
        super().__init__()
        self.close()
        self.setGeometry(20, 500, 500, 400)
        self.setWindowTitle("Response Surface, Linear Model with Interactions")
        self.setFixedSize(500, 400)
        plot_first_order_canvas(self)
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setFixedSize(500, 20)
        self.addToolBar(bar)
        save_fig = QAction(QIcon("/PATH/disk-black.png"), "Save Figure", self)
        save_fig.setCheckable(True)
        save_fig.triggered.connect(self.save_surface_1_figure)
        bar.addAction(save_fig)
        self.show()
    # Defining Save Figure Surface 1
    def save_surface_1_figure(self):
        surface_1_path = QFileDialog.getSaveFileName(self, "Save Figure", "*.jpg")
        file_name = str(surface_1_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            plt.figure(2).savefig(file_name_1)
        return
# First Order Projection Plot Canvas
class plot_first_order_projection_canvas(FigureCanvas):
    # Defining the buffer and variables in the buffer as globals
    global temp_input_data, f_1, f_2, resp_f, const, results_buffer_1
    def __init__(self, parent):
        processing()
        first_order_fitting()
        fact_names()
        plt.close(3)
        my_results = pd.read_csv(io.StringIO(results_buffer_1.getvalue()))
        const = my_results[" Coefficient  "].explode()
        s1 = (max(f_1) - min(f_1)) / 100
        s2 = (max(f_2) - min(f_2)) / 100
        x = np.arange(min(f_1), max(f_1), s1)
        y = np.arange(min(f_2), max(f_2), s2)
        z = const[0] + const[1] * x + const[2] * y + const[3] * x * y
        X, Y = np.meshgrid(x, y)
        Z = const[0] + const[1] * X + const[2] * Y + const[3] * X * Y
        d = np.array([x, y, z])
        d = np.transpose(d)
        planes = []
        for i in range(0, len(z)):
            planes.append(float(d[i, 2]))
        planes = sorted(planes)
        self.sc_plot = plt.figure(3, figsize=(5, 4))
        self.eje1 = plt.axes()
        self.eje1.set_title("Surface Projection, Linear Model with Interactions")
        self.eje1.set_xlabel(f_1_name)
        self.eje1.set_ylabel(f_2_name)
        plt.contour(X, Y, Z, levels=planes[0:-1], cmap="winter")
        plt.contour(X, Y, Z, colors="black")
        plt.draw()
        super().__init__(self.sc_plot)
        self.setParent(parent)
# First Order Projection Plot Main
class plot_first_order_projection(QMainWindow):
    def __init__(self):
        super().__init__()
        self.close()
        self.setGeometry(700, 500, 500, 400)
        self.setWindowTitle("Surface Projection, Linear Model with Interactions")
        self.setFixedSize(500, 400)
        plot_first_order_projection_canvas(self)
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setFixedSize(500, 20)
        self.addToolBar(bar)
        save_fig = QAction(QIcon("/PATH/disk-black.png"), "Save Figure", self)
        save_fig.setCheckable(True)
        save_fig.triggered.connect(self.save_surface_1_projection)
        bar.addAction(save_fig)
        self.show()
# Defining Save Figure Surface 1
    def save_surface_1_projection(self):
        surface_1_proj = QFileDialog.getSaveFileName(self, "Save Figure", "*.jpg")
        file_name = str(surface_1_proj)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            plt.figure(3).savefig(file_name_1)
        return

# Quadratic Plot Canvas
class plot_quadratic_canvas(FigureCanvas):
    # Defining the buffer and variables in the buffer as globals
    global temp_input_data, f_1, f_2, resp_f, const, results_buffer
    def __init__(self, parent):
        processing()
        second_order_fitting()
        fact_names()
        plt.close(4)
        my_results = pd.read_csv(io.StringIO(results_buffer.getvalue()))
        const = my_results[" Coefficient  "].explode()
        self.sc_plot = plt.figure(4, figsize=(5, 4))
        self.eje1 = plt.axes(projection="3d")
        self. eje1.scatter(f_1, f_2, resp_f, s=50, c="red")
        self.eje1.set_title("Response Surface, Quadratic Model with Interactions")
        self.eje1.set_xlabel(f_1_name)
        self.eje1.set_ylabel(f_2_name)
        self.eje1.set_zlabel(resp_name)
        s1 = (max(f_1) - min(f_1)) / 100
        s2 = (max(f_2) - min(f_2)) / 100
        x = np.arange(min(f_1), max(f_1), s1)
        y = np.arange(min(f_2), max(f_2), s2)
        X, Y = np.meshgrid(x, y)
        z = const[0] + const[1] * x + const[2] * y + const[3] * x ** 2 + const[
            4] * y ** 2 + const[5] * x * y + const[6] * x ** 2 * y ** 2
        Z = const[0] + const[1] * X + const[2] * Y + const[3] * X ** 2 + const[
            4] * Y ** 2 + const[5] * X * Y + const[6] * X ** 2 * Y ** 2
        self.eje1.plot_surface(X, Y, Z, cmap="winter")
        plt.draw()
        super().__init__(self.sc_plot)
        self.setParent(parent)

# Quadratic Plot Main
class plot_quadratic(QMainWindow):
    def __init__(self):
        super().__init__()
        self.close()
        self.setGeometry(1150, 20, 500, 400)
        self.setWindowTitle("Response Surface, Quadratic Model with Interactions")
        self.setFixedSize(500, 400)
        plot_quadratic_canvas(self)
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setFixedSize(500, 20)
        self.addToolBar(bar)
        save_fig = QAction(QIcon("/PATH/disk-black.png"), "Save Figure", self)
        save_fig.setCheckable(True)
        save_fig.triggered.connect(self.save_surface_2_figure)
        bar.addAction(save_fig)
        self.show()
# Defining Save Figure Surface 2
    def save_surface_2_figure(self):
        surface_2_path = QFileDialog.getSaveFileName(self, "Save Figure", "*.jpg")
        file_name = str(surface_2_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            plt.figure(4).savefig(file_name_1)
        return

# Quadratic Projection Plot Canvas
class plot_quadratic_projection_canvas(FigureCanvas):
    # Defining the buffer and variables in the buffer as globals
    global temp_input_data, f_1, f_2, resp_f, const, results_buffer
    def __init__(self, parent):
        processing()
        second_order_fitting()
        fact_names()
        plt.close(5)
        my_results = pd.read_csv(io.StringIO(results_buffer.getvalue()))
        const = my_results[" Coefficient  "].explode()
        s1 = (max(f_1) - min(f_1)) / 100
        s2 = (max(f_2) - min(f_2)) / 100
        x = np.arange(min(f_1), max(f_1), s1)
        y = np.arange(min(f_2), max(f_2), s2)
        X, Y = np.meshgrid(x, y)
        z = const[0] + const[1] * x + const[2] * y + const[3] * x ** 2 + const[
            4] * y ** 2 + const[5] * x * y + const[6] * x ** 2 * y ** 2
        Z = const[0] + const[1] * X + const[2] * Y + const[3] * X ** 2 + const[
            4] * Y ** 2 + const[5] * X * Y + const[6] * X ** 2 * Y ** 2
        d = np.array([x, y, z])
        d = np.transpose(d)
        planes = []
        for i in range(0, len(z)):
            planes.append(float(d[i, 2]))
        planes = sorted(planes)
        self.sc_plot = plt.figure(5, figsize=(5, 4))
        self.eje1 = plt.axes()
        self.eje1.set_title("Surface Projection, Quadratic Model with Interactions")
        self.eje1.set_xlabel(f_1_name)
        self.eje1.set_ylabel(f_2_name)
        plt.contour(X, Y, Z, levels=planes[0:-1], cmap="winter")
        plt.contour(X, Y, Z, colors="black")
        plt.draw()
        super().__init__(self.sc_plot)
        self.setParent(parent)
# First Order Projection Plot Main
class plot_quadratic_projection(QMainWindow):
    def __init__(self):
        super().__init__()
        self.close()
        self.setGeometry(1150, 500, 500, 400)
        self.setWindowTitle("Surface Projection, Quadratic Model with Interactions")
        self.setFixedSize(500, 400)
        plot_quadratic_projection_canvas(self)
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setFixedSize(500, 20)
        self.addToolBar(bar)
        save_fig = QAction(QIcon("/PATH/disk-black.png"), "Save Figure", self)
        save_fig.setCheckable(True)
        save_fig.triggered.connect(self.save_surface_2_projection)
        bar.addAction(save_fig)
        self.show()
# Defining Save Figure Surface 2 Projection
    def save_surface_2_projection(self):
        surface_2_projection = QFileDialog.getSaveFileName(self, "Save Figure", "*.jpg")
        file_name = str(surface_2_projection)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            plt.figure(5).savefig(file_name_1)
        return
#-----------------------------------------------------------------------------------------------------------------------
# Linear Model with interactions
class first_order_fitting(QTableWidget):
    def __init__(self):
        global det_coeff, Fisher_Val, Fisher_prob, results_buffer_1, f_1_name, f_2_name, resp_name
        # Init
        super().__init__()
        self.calculations_first_order()
        fact_names()
        font = QFont()
        font.setPointSize(15)
        font.setBold(True)
        # Create button to export table
        save_data = QPushButton(self)
        save_data.setIcon(QIcon("/PATH/disk-black.png"))
        save_data.setCheckable(True)
        save_data.setAccessibleDescription("Export as CSV")
        save_data.clicked.connect(self.save_table_1)
        self.setFont(font)
        self.setWindowTitle("Fitting Data to linear Model with interactions")
        self.setGeometry(800, 100, 900, 300)
        self.setFixedSize(980, 300)
        self.setColumnCount(5)
        self.setColumnWidth(0, 110)
        self.setColumnWidth(1, 110)
        self.setColumnWidth(2, 150)
        self.setColumnWidth(3, 220)
        self.setColumnWidth(4, 220)
        self.setRowCount(5)
        self.setAlternatingRowColors(True)
        c_header = ["Coefficient", "t-Value", "Probability (t)", "Conf Int(95%), lower limit", "Conf Int(95%), higher limit"]
        f_header = ["Parameter", "Constant", f_1_name, f_2_name, f_1_name +"*"+ f_2_name]
        c_header_i = []
        f_header_i = []
        for i in range(0, len(c_header)):
            c_header_i.append(QTableWidgetItem(c_header[i]))
        for j in range(0, len(f_header)):
            f_header_i.append(QTableWidgetItem(f_header[j]))
        for x in range(0, len(c_header_i)):
            self.setHorizontalHeaderItem(x, c_header_i[x])
            c_header_i[x].setTextAlignment(Qt.AlignCenter)
        for y in range(0, len(f_header_i)):
            self.setVerticalHeaderItem(y, f_header_i[y])
            f_header_i[y].setTextAlignment(Qt.AlignCenter)
        self.r_squared = QLabel(self)
        self.f_value = QLabel(self)
        self.prob_f_value = QLabel(self)
        self.r_squared.setGeometry(180, 150, 180, 100)
        self.r_squared.setText("R-squared = " + str(round(det_coeff, 4)))
        self.f_value.setGeometry(180, 180, 180, 100)
        self.f_value.setText("F-value = " + str(round(Fisher_Val, 4)))
        self.prob_f_value.setGeometry(180, 210, 180, 100)
        self.prob_f_value.setText("Probability(F) = " + str(round(Fisher_prob, 10)))
        # Extracting data from file to show in screen
        my_results = pd.read_csv(io.StringIO(results_buffer_1.getvalue()))
        # Probability associated to t value
        t_prob_item_1 = my_results["  Probability (t)  "].explode()
        tp_item_1 = []
        # t- value
        t_val_1 = my_results["  t value  "].explode()
        t_item_1 = []
        # Coefficients
        coef_val_1 = my_results[" Coefficient  "].explode()
        coef_item_1 = []
        # Limits of CI
        lower_CI = my_results["low_CI"].explode()
        lower_CI_item = []
        higher_CI = my_results["high_CI"].explode()
        higher_CI_item = []
        # Loop to show data in table
        for u in range(0, len(t_prob_item_1)):
            # t_probability
            tp_item_1.append(QTableWidgetItem(str(round(t_prob_item_1[u], 10))))
            # t-value
            t_item_1.append(QTableWidgetItem(str(round(t_val_1[u], 3))))
            # Coefficients
            coef_item_1.append(QTableWidgetItem(str(round(coef_val_1[u], 3))))
            # Limits of CI
            lower_CI_item.append(QTableWidgetItem(str(round((lower_CI[u]), 3))))
            higher_CI_item.append(QTableWidgetItem(str(round((higher_CI[u]), 3))))
        # Loop for inserting items in Table
        for x in range(0, len(tp_item_1)):
            # Setting the item flags to static editable values in table
            tp_item_1[x].setFlags(tp_item_1[x].flags() ^ Qt.ItemIsEditable)
            t_item_1[x].setFlags(t_item_1[x].flags() ^ Qt.ItemIsEditable)
            coef_item_1[x].setFlags(coef_item_1[x].flags() ^ Qt.ItemIsEditable)
            lower_CI_item[x].setFlags(lower_CI_item[x].flags() ^ Qt.ItemIsEditable)
            higher_CI_item[x].setFlags(higher_CI_item[x].flags() ^ Qt.ItemIsEditable)
            self.setItem(x + 1, 2, tp_item_1[x])
            self.setItem(x + 1, 1, t_item_1[x])
            self.setItem(x + 1, 0, coef_item_1[x])
            self.setItem(x + 1, 3, lower_CI_item[x])
            self.setItem(x + 1, 4, higher_CI_item[x])
            # Centering Items
            tp_item_1[x].setTextAlignment(Qt.AlignCenter)
            t_item_1[x].setTextAlignment(Qt.AlignCenter)
            coef_item_1[x].setTextAlignment(Qt.AlignCenter)
            higher_CI_item[x].setTextAlignment(Qt.AlignCenter)
            lower_CI_item[x].setTextAlignment(Qt.AlignCenter)
    # ----------------------------------------------------------------------------------------------------------
    # Regression First Order
    def calculations_first_order(self):
        global f_1, f_2, resp_f, results_buffer_1, det_coeff, Fisher_Val, Fisher_prob, Results, results_buffer_1, f_1_name, f_2_name, resp_name
        # Extracting columns data from Data Frame stored in buffer (temp_input_data)
        extracting_data_from_data_frame()
        # Obtaining other terms of the model (interaction, and unit column)
        c = []
        inter_f1_f2 = []
        # Loop to Matrix
        for h in range(0, len(f_1)):
            c.append(1)
            inter_f1_f2.append(f_1[h] * f_2[h])
        # Matrix Creation
        M = np.array([c, f_1, f_2, inter_f1_f2])
        M = np.transpose(M)
        # Regression
        # Modelling process
        model = sm.OLS(resp_f, M)
        regression = model.fit()
        # Table with output results
        det_coeff = regression.rsquared
        Fisher_Val = regression.fvalue
        Fisher_prob = regression.f_pvalue
        t_value = []
        t_prob = []
        constants = []
        ic_0 = []
        ic_1 = []
        interval_c = regression.conf_int()
        for i in range(0, len(interval_c)):
            ic_0.append(interval_c[0][i])
            ic_1.append(interval_c[1][i])
        for k in range(0, 4):
            constants.append(round((regression.params[k]), 4))
            t_value.append(round((regression.tvalues[k]), 3))
            t_prob.append(round((regression.pvalues[k]), 3))
        Results = {"  Parameter  ": ["Constant", f_1_name, f_2_name, f_1_name + "*" + f_2_name],
                   " Coefficient  ": constants, "  t value  ": t_value,
                   "  Probability (t)  ": t_prob, "low_CI": ic_0, "high_CI": ic_1,
                   " R-squared": [det_coeff, "", "", ""],
                   "F-Value": [Fisher_Val, "", "", ""],
                   "Probability(F)": [Fisher_prob, "", "", ""]}
        Results = pd.DataFrame(Results, index=None)
        results_buffer_1 = io.StringIO()
        Results.to_csv(results_buffer_1)
    # Save_table
    def save_table_1(self):
        global Results
        file_path = QFileDialog.getSaveFileName(self, "Save File ", "*.csv")
        file_name = str(file_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            Results.to_csv(file_name_1, index=None)
        return
# ----------------------------------------------------------------------------------------------------------------------
# Class Fitting to Quadratic Model
class second_order_fitting(QTableWidget):
    def __init__(self):
        global det_coeff, Fisher_Val, Fisher_prob, results_buffer_1, coef_val_2, my_results, max_value, f_1max, f_2max
        # Init
        super().__init__()
        self.close()
        self.calculations_second_order()
        font = QFont()
        font.setPointSize(15)
        font.setBold(True)
        # Create button to export table
        save_data = QPushButton(self)
        save_data.setIcon(QIcon("/PATH/disk-black.png"))
        save_data.setCheckable(True)
        save_data.setAccessibleDescription("Export as CSV")
        save_data.clicked.connect(self.save_table_2)
        self.setFont(font)
        self.setWindowTitle("Fitting Data to quadratic model with interactions")
        self.setGeometry(800, 500, 1000, 400)
        self.setFixedSize(1000, 400)
        self.setColumnCount(5)
        self.setColumnWidth(0, 110)
        self.setColumnWidth(1, 110)
        self.setColumnWidth(2, 150)
        self.setColumnWidth(3, 220)
        self.setColumnWidth(4, 220)
        self.setRowCount(8)
        self.setAlternatingRowColors(True)
        c_header_2 = ["Coefficient", "t-Value", "Probability (t)", "Conf Int(95%), lower limit", "Conf Int(95%), higher limit"]
        f_header_2 = ["Parameter", "Constant", f_1_name, f_2_name, f_1_name+"^2", f_2_name+"^2",f_1_name+"*"+f_2_name,
                      "("+f_1_name+"*"+f_2_name+")^2"]
        c_header_i_2 = []
        f_header_i_2 = []
        for i in range(0, len(c_header_2)):
            c_header_i_2.append(QTableWidgetItem(c_header_2[i]))
        for j in range(0, len(f_header_2)):
            f_header_i_2.append(QTableWidgetItem(f_header_2[j]))
        for x in range(0, len(c_header_i_2)):
            self.setHorizontalHeaderItem(x, c_header_i_2[x])
            c_header_i_2[x].setTextAlignment(Qt.AlignCenter)
        for y in range(0, len(f_header_i_2)):
            self.setVerticalHeaderItem(y, f_header_i_2[y])
            f_header_i_2[y].setTextAlignment(Qt.AlignCenter)
        self.r_squared = QLabel(self)
        self.f_value = QLabel(self)
        self.prob_f_value = QLabel(self)
        self.f1_max = QLabel(self)
        self.f2_max = QLabel(self)
        self.max_resp = QLabel(self)
        self.r_squared.setGeometry(220, 250, 180, 100)
        self.r_squared.setText("R-squared = " + str(round(det_coeff, 4)))
        self.f_value.setGeometry(220, 280, 180, 100)
        self.f_value.setText("F-value = " + str(round(Fisher_Val, 4)))
        self.prob_f_value.setGeometry(220, 310, 180, 100)
        self.prob_f_value.setText("Probability(F) = " + str(round(Fisher_prob, 10)))
        self.max_resp.setGeometry(630, 250, 220, 100)
        self.max_resp.setText("Maximum Value = " + str(round(max_value, 2)))
        self.f1_max.setGeometry(630, 280, 220, 100)
        self.f1_max.setText(f_1_name + " = " + str(round(f_1max, 2)))
        self.f2_max.setGeometry(630, 310, 220, 100)
        self.f2_max.setText(f_2_name + " = " + str(round(f_2max, 2)))
        # Extracting data from file to show in screen--> check how to do it with a memory buffer
        my_results = pd.read_csv(io.StringIO(results_buffer.getvalue()))
        # Probability associated to t value
        t_prob_item_2 = my_results["  Probability (t)  "].explode()
        tp_item_2 = []
        # t- value
        t_val_2 = my_results["  t value  "].explode()
        t_item_2 = []
        # Coefficients
        coef_val_2 = my_results[" Coefficient  "].explode()
        coef_item_2 = []
        # Limits of CI
        lower_CI_2 = my_results["low_CI"].explode()
        lower_CI_item_2 = []
        higher_CI_2 = my_results["high_CI"].explode()
        higher_CI_item_2 = []
        # Loop to show data in table
        for u in range(0, len(t_prob_item_2)):
            # t_probability
            tp_item_2.append(QTableWidgetItem(str(round(t_prob_item_2[u], 10))))
            # t-value
            t_item_2.append(QTableWidgetItem(str(round(t_val_2[u], 4))))
            # Coefficients
            coef_item_2.append(QTableWidgetItem(str(round(coef_val_2[u], 4))))
            # Limits of CI
            lower_CI_item_2.append(QTableWidgetItem(str(round((lower_CI_2[u]), 3))))
            higher_CI_item_2.append(QTableWidgetItem(str(round((higher_CI_2[u]), 3))))
        # Loop for inserting items in Table
        for x in range(0, len(tp_item_2)):
            # Setting the item flags to static editable values in table
            tp_item_2[x].setFlags(tp_item_2[x].flags() ^ Qt.ItemIsEditable)
            t_item_2[x].setFlags(t_item_2[x].flags() ^ Qt.ItemIsEditable)
            coef_item_2[x].setFlags(coef_item_2[x].flags() ^ Qt.ItemIsEditable)
            lower_CI_item_2[x].setFlags(lower_CI_item_2[x].flags() ^ Qt.ItemIsEditable)
            higher_CI_item_2[x].setFlags(higher_CI_item_2[x].flags() ^ Qt.ItemIsEditable)
            # Inserting items in table
            self.setItem(x + 1, 2, tp_item_2[x])
            self.setItem(x + 1, 1, t_item_2[x])
            self.setItem(x + 1, 0, coef_item_2[x])
            self.setItem(x + 1, 3, lower_CI_item_2[x])
            self.setItem(x + 1, 4, higher_CI_item_2[x])
            # Centering Items
            tp_item_2[x].setTextAlignment(Qt.AlignCenter)
            t_item_2[x].setTextAlignment(Qt.AlignCenter)
            coef_item_2[x].setTextAlignment(Qt.AlignCenter)
            higher_CI_item_2[x].setTextAlignment(Qt.AlignCenter)
            lower_CI_item_2[x].setTextAlignment(Qt.AlignCenter)

    # Regression Second Order
    def calculations_second_order(self):
        global f_1, f_2, resp_f, results_buffer, det_coeff, Fisher_Val, Fisher_prob, max_value, f_1max, f_2max, Results_2, f_1_max, f_2_max
        # Extracting column data from Data Frame stored in buffer (temp_input_data)
        extracting_data_from_data_frame()
        # Obtaining other terms of the model (interaction, quadratic terms and unit column)
        c = []
        inter_f1_f2 = []
        factor_1_quad = []
        factor_2_quad = []
        inter_f1_f_2_quad = []
        # Loop to Matrix
        for h in range(0, len(f_1)):
            c.append(1)
            inter_f1_f2.append(f_1[h] * f_2[h])
            factor_1_quad.append((f_1[h]) ** 2)
            factor_2_quad.append((f_2[h]) ** 2)
            inter_f1_f_2_quad.append((f_1[h] ** 2) * (f_2[h] ** 2))
        # Matrix Creation
        M = np.array([c, f_1, f_2, factor_1_quad, factor_2_quad, inter_f1_f2, inter_f1_f_2_quad])
        M = np.transpose(M)
        # Regression
        # Modelling process
        model = sm.OLS(resp_f, M)
        regression = model.fit()
        # Table with output results
        det_coeff = regression.rsquared
        Fisher_Val = regression.fvalue
        Fisher_prob = regression.f_pvalue
        t_value = []
        t_prob = []
        constants = []
        ic_0_2 = []
        ic_1_2 = []
        interval_c_2 = regression.conf_int()
        for i in range(0, len(interval_c_2)):
            ic_0_2.append(interval_c_2[0][i])
            ic_1_2.append(interval_c_2[1][i])
        for k in range(0, 7):
            constants.append(regression.params[k])
            t_value.append(regression.tvalues[k])
            t_prob.append(regression.pvalues[k])
        Results_2 = {"  Factor  ": ["Constant", "Factor_1", "Factor_2", "Factor_1^2", "Factor_2^2", "(Factor_1*Factor_2)",
                                  "(Factor_1*Factor_2)^2"],
                   " Coefficient  ": constants,
                   "  t value  ": t_value,
                   "  Probability (t)  ": t_prob, "low_CI": ic_0_2, "high_CI": ic_1_2,
                    " R-squared": [det_coeff, "", "", "", "", "", ""],
                     "F-Value": [Fisher_Val, "", "", "", "", "", ""],
                     "Probability(F)": [Fisher_prob, "", "", "", "", "", ""]
                     }
        Results_2 = pd.DataFrame(Results_2, index=None)
        results_buffer = io.StringIO()
        Results_2.to_csv(results_buffer)
        # Optimization * this is a rough analysis. Need to implement reliable numerical methods to optimize
        s1 = (max(f_1) - min(f_1)) / 50
        s2 = (max(f_2) - min(f_2)) / 50
        x = np.arange(min(f_1), max(f_1), s1)
        y = np.arange(min(f_2), max(f_2), s2)
        X = []; Y = []; Z= []; e = []; f = []
        for i in range(0, 50):
            X.append(x[i])
            Y.append(y[i])
        for a in range(0, len(X)):
            for b in range(0, len(Y)):
                d = np.array([X[a], Y[b]])
                d = np.transpose(d)
                e.append(d[0])
                f.append(d[1])
        for h in range(0, len(e)):
            Z.append(constants[0] + constants[1] * e[h] + constants[2] * f[h] + constants[3] * e[h] ** 2
                     + constants[4] * f[h] ** 2 + constants[5] * e[h] * f[h] + constants[6] * e[h] ** 2 * f[h] ** 2)
        for p in range(0, len(Z)):
            if Z[p] == max(Z):
                f_1max = e[p]
                f_2max = f[p]
                max_value = Z[p]
    # Save_table
    def save_table_2(self):
        global Results_2
        file_path = QFileDialog.getSaveFileName(self, "Save File ", "*.csv")
        file_name = str(file_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            Results_2.to_csv(file_name_1, index=None)
        return

#================================================================================================================
# Input Table
class input_table(QTableWidget):
    global temp_input_data, f_1_name, f_2_name, resp_name, header_change_0, header_change_1, header_change_2
    def __init__(self, r, c):
        global temp_input_data, f_1_name, f_2_name, resp_name, header_change_0, header_change_1, header_change_2
        super().__init__(100, 3)
        self.iw = None
        self.bm = None
        self.fn = None
        self.w = None
        self.w1 = None
        self.ps = None
        self.fo = None
        self.fop = None
        self.qo = None
        self.qop = None
        self.ne = None
        self.clear_all()
        font = QFont()
        font.setPointSize(12)
        self.setFont(font)
        self.setWindowTitle("Bifactorial Design ver. 1.0")
        self.setGeometry(800, 100, 680, 650)
        self.setFixedSize(720, 500)
        col_headers = ["FACTOR 1", 'FACTOR 2', 'RESPONSE']
        self.setHorizontalHeaderLabels(col_headers)
        self.setCurrentCell(0, 0)
        self.init_ui()
        # Presentation
        present = QLabel(self)
        present.setGeometry(400, 40, 350, 15)
        present.setText("BIFACTORIAL DESIGN ver. 1.0 (2022)")
        present_1 = QLabel(self)
        present_1.setGeometry(380, 60, 350, 15)
        present_1.setText("An open source project by Armando Hernández")
        # Logo
        logo_im = QLabel(self)
        logo_im.setGeometry(330, 80, 390, 500)
        logo_im.setAutoFillBackground(True)
        logo_im.setAlignment(Qt.AlignCenter)
        logo_im.setPixmap(QPixmap("/PATH/Presentation_test_2.png"))
    # Def Init_UI
    def init_ui(self):
        self.cellChanged.connect(self.c_current)
        self.cellEntered.connect(self.c_current)
        self.cellClicked.connect(self.c_current)
        self.cellPressed.connect(self.c_current)
        self.cellDoubleClicked.connect(self.c_current)
        # Context Menu Options
        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        #quit_option = QAction("Quit", self)
        select_option = QAction("Select", self)
        copy_option = QAction("Copy", self)
        paste_option = QAction("Paste", self)
        #quitAction.triggered.connect(qApp.quit())
        self.addActions([select_option, copy_option, paste_option])
    # Def Current Cell
    def c_current(self):
        global temp_input_data, result, fact_1, test
        result = []
        for col in range(0, self.columnCount()):
            rows = []
            for row in range(0, self.rowCount()):
                value = self.item(row, col)
                if value:
                    rows.append(value.text())
                    value.setTextAlignment(Qt.AlignCenter)
                    if value.text() == "":
                        rows.clear()
            result.append(rows)
        fact_1 = result[0]
        fact_2 = result[1]
        resp = result[2]
        test = np.array([fact_1, fact_2, resp], dtype=object)
        test = np.transpose(test)
        test_data = pd.DataFrame(test)
        # Saving to buffer
        temp_input_data = io.StringIO()
        test_data.to_csv(temp_input_data)
    # ================================================================================================================
    def data_in_table_for_plot_scattered(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_scattered_connection()
            return
    # Def Data in table for linear Model with interactions
    def data_in_table_for_fit_first_order(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.first_order_fitting_c()
            return

    # Def Data in table for quadratic Model with interactions
    def data_in_table_for_fit_quadratic(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.second_order_fitting_c()
            return

    # Def Data in table for surface linear Model
    def data_in_table_for_surface_first_order(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_first_order_connection()
            return

    # Def Data in table for surface quadratic Model
    # Plot Surface Decision if no- data entry, then pop-up something...e.g. "You haven't entered data yet..."
    def data_in_table_for_surface_quadratic(self):
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            self.plot_quadratic_connection()
            return

    # Window to connect the Name of Factors
    def factors_names_connection(self):
        global f_1_name, f_2_name, resp_name, header_change_0, header_change_1, header_change_2, f_1_name_temp, \
            f_2_name_temp, resp_name_temp
        # Change of Name Factor 1
        f_1_name_temp = f_1_name
        f_1_name, header_change_0 = QInputDialog.getText(self, "Factor 1 Name", 'FACTOR 1:', QLineEdit.Normal, f_1_name)
        if f_1_name == "":
            f_1_name = f_1_name_temp
        if header_change_0:
            self.setHorizontalHeaderItem(0, QTableWidgetItem(f_1_name))
        # Change of Name Factor 2
        f_2_name_temp = f_2_name
        f_2_name, header_change_1 = QInputDialog.getText(self, "Factor 2 Name", 'FACTOR 2:', QLineEdit.Normal, f_2_name)
        if f_2_name == "":
            f_2_name = f_2_name_temp
        if header_change_1:
            self.setHorizontalHeaderItem(1, QTableWidgetItem(f_2_name))
        # Change of Name Factor 3
        resp_name_temp = resp_name
        resp_name, header_change_2 = QInputDialog.getText(self, "Response Name", "RESPONSE:", QLineEdit.Normal,
                                                          resp_name)
        if resp_name == "":
            resp_name = resp_name_temp
        if header_change_2:
            self.setHorizontalHeaderItem(2, QTableWidgetItem(resp_name))
        self.c_current()
        return

        # Window to pop up "Not enough data"
    def not_enough_data_connection(self):
        global temp_input_data
        if self.ne is None:
            self.ne = not_enough_data()
            self.ne.show()
        elif not_enough_data():
            self.ne = not_enough_data()
            self.ne.show()
        else:
            self.ne.close()
            self.ne = None

    # Window to show the Scatter Plot
    def plot_scattered_connection(self):
        global temp_input_data
        if self.ps is None:
            self.ps = plot_scattered()
            self.ps.show()
        elif plot_scattered():
            self.ps = plot_scattered()
            self.ps.show()
        else:
            self.ps.close()
            self.ps = None
    # Window to show the First -Order Plot
    def plot_first_order_connection(self):
        if self.fo is None and self.fop is None:
            self.fo = plot_first_order()
            self.fop = plot_first_order_projection()
            self.fo.show()
            self.fop.show()
        elif plot_first_order():
            self.fo = plot_first_order()
            self.fop = plot_first_order_projection()
            self.fo.show()
            self.fop.show()
        else:
            self.fo.close()
            self.fo = None
            self.fop.close()
            self.fop = None

    # Window to show the Second-Order Plot
    def plot_quadratic_connection(self):
        if self.qo is None and self.qop is None:
            self.qo = plot_quadratic()
            self.qop = plot_quadratic_projection()
            self.qo.show()
            self.qop.show()
        elif plot_first_order():
            self.qo = plot_quadratic()
            self.qop = plot_quadratic_projection()
            self.qo.show()
            self.qop.show()
        else:
            self.qo.close()
            self.qop = None
            self.qo.close()
            self.qop = None
    # Window to show the regression results regarding linear fitting
    def first_order_fitting_c(self):
        if self.w is None:
            self.w = first_order_fitting()
            self.w.show()
        elif first_order_fitting():
            self.w = first_order_fitting()
            self.w.show()
        else:
            self.w.close()
            self.w = None

    # Window to show the regression results regarding quadratic fitting
    def second_order_fitting_c(self):
        if self.w1 is None:
            self.w1 = second_order_fitting()
            self.w1.show()
        elif second_order_fitting():
            self.w1 = second_order_fitting()
            self.w1.show()
        else:
            self.w1.close()
            self.w1 = None

    # Open file
    def open_file(self):
        global temp_input_data, f_1_name, f_2_name, resp_name
        uploaded_data = ""
        file_path = QFileDialog.getOpenFileName(self, "Open File ", "d:\\", "*.csv *.xls *.xlsx")
        file_name = str(file_path)
        str_list = []
        c = 1
        file_name_1 = ""
        for i in range(1, len(file_name)):
            str_list.append(file_name[i])
        for j in range(0, len(str_list)):
            while str_list[c] != "'":
                file_name_1 = file_name_1 + str_list[c]
                c = c + 1
        if file_name_1 != "":
            if file_name_1[-1] == "v":
                uploaded_data = pd.read_csv(file_name_1)
            elif file_name_1[-1] == "s" or file_name_1[-1] == "x":
                uploaded_data = pd.read_excel(file_name_1)
            f_1_name = uploaded_data.columns[0]
            f_2_name = uploaded_data.columns[1]
            resp_name = uploaded_data.columns[2]
            # Removing Headers
            modify_uploaded_data = uploaded_data.rename(
                columns={uploaded_data.columns[0]: "0", uploaded_data.columns[1]: "1", uploaded_data.columns[2]: "2"})
            # Clearing the Table and the previous results
            self.clear_cells_input_table()
            # Extracting Data
            c_0 = modify_uploaded_data["0"].explode()
            c_1 = modify_uploaded_data["1"].explode()
            c_2 = modify_uploaded_data["2"].explode()
            c_0_item = []
            c_1_item = []
            c_2_item = []
            for u in range(0, len(c_2)):
                c_0_item.append(QTableWidgetItem(str(c_0[u])))
                c_1_item.append(QTableWidgetItem(str(c_1[u])))
                c_2_item.append(QTableWidgetItem(str(c_2[u])))
            for x in range(len(c_2_item)):  # Setting Data in Table
                self.setItem(x, 0, c_0_item[x])
                self.setItem(x, 1, c_1_item[x])
                self.setItem(x, 2, c_2_item[x])
            self.setHorizontalHeaderLabels([f_1_name, f_2_name, resp_name])
            return
        else:
            return

    # Saving Data as CSV
    def save_data(self):
        global temp_input_data, result, f_1_name, f_2_name, resp_name
        self.c_current()
        if len(test) <= 3:
            self.not_enough_data_connection()
            return
        else:
            # Temporary reading data from buffer to save it later
            temp_read = pd.read_csv(io.StringIO(temp_input_data.getvalue()))
            file_path = QFileDialog.getSaveFileName(self, "Save File ", "*.csv")
            file_name = str(file_path)
            str_list = []
            c = 1
            file_name_1 = ""
            for i in range(1, len(file_name)):
                str_list.append(file_name[i])
            for j in range(0, len(str_list)):
                while str_list[c] != "'":
                    file_name_1 = file_name_1 + str_list[c]
                    c = c + 1
            if file_name_1 != "":
                modify_temp_read = temp_read.rename(
                    columns={temp_read.columns[1]: f_1_name, temp_read.columns[2]: f_2_name,
                             temp_read.columns[3]: resp_name})  # Replacing the headers names
                modify_temp_read.to_csv(file_name_1, index=None, columns=[f_1_name, f_2_name, resp_name])
            return

    # Refresh Cells
    def clear_cells_input_table(self):
        global temp_input_data, f_1_name, f_2_name, resp_name, resp_name_temp, f_1_name_temp, f_2_name_temp
        self.clear()
        fact_names()
        col_headers = [f_1_name, f_2_name, resp_name]
        self.setHorizontalHeaderLabels(col_headers)
        temp_input_data = io.StringIO()

    # Clear All
    def clear_all(self):
        global temp_input_data, f_1_name, f_2_name, resp_name
        # closing scatter-plot 3D-figure
        if self.ps:
            self.ps.close()
            self.ps = None
        # closing First-order Plot 3D-figure
        if self.fo:
            self.fo.close()
            self.fo = None
        if self.fop:
            self.fop.close()
            self.fop = None
        # closing Quadratic Plot 3D-figure
        if self.qo:
            self.qo.close()
            self.qo = None
        if self.qop:
            self.qop.close()
            self.qop = None
        # closing first-order results
        if self.w is not None:
            self.w.clear()
            self.w.close()
            self.w = None
        # closing quadratic fitting results
        if self.w1 is not None:
            self.w1.clear()
            self.w1.close()
            self.w1 = None
        temp_input_data = io.StringIO()
        self.clear()
        col_headers = ["FACTOR 1", "FACTOR 2", "RESPONSE"]
        self.setHorizontalHeaderLabels(col_headers)

# Main Script


if __name__ == "__main__":
    factorial_exp = QApplication(sys.argv)
    fact_design = initial_window()
    fact_design.show()
    sys.exit(factorial_exp.exec_())