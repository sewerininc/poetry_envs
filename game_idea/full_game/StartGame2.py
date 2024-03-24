import sys
import random
import time
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QLineEdit, QLabel, QMessageBox, QGraphicsView, QProgressBar
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
from IPython import display
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import tempfile
import os
import numpy as np
from matplotlib.ticker import MaxNLocator


# Write everything into the class

class GameApp(QWidget):
    def __init__(self):
        super().__init__()

        self.names = []
        self.name_count = {}
        self.initUI()

    def initUI(self):
        # testing
        # self.names = ["Alice", "Bob", "Charlie", "David", "Eve"]
        # self.name_count = {name: 0 for name in self.names}
        
        self.setWindowTitle('Drinking Game')
        self.setGeometry(100, 100, 400, 200)

        self.name_label = QLabel('Enter a names: (format: mikkel tobais julie tomas)')
        self.name_input = QLineEdit()
        self.add_names_button = QPushButton('Add')
        self.run_button = QPushButton('Run')
        self.progress_bar = QProgressBar()
        
        self.add_names_button.clicked.connect(self.add_names)
        self.run_button.clicked.connect(self.run_game)
        self.figure_canvas = FigureCanvasQTAgg(plt.figure())

        vbox = QVBoxLayout()
        vbox.addWidget(self.name_label)
        vbox.addWidget(self.name_input)
        vbox.addWidget(self.add_names_button)
        vbox.addWidget(self.run_button)
        vbox.addWidget(self.progress_bar)
        vbox.addWidget(self.figure_canvas)

        self.setLayout(vbox)
    
    
    def add_names(self):
        try:
            all_names = self.name_input.text().split(" ")
            if all_names == [""]:
                self.show_error_message("Wrong Format for lists no names added")
                
                return
            
            for name in all_names:
                self.names.append(name)

            self.name_count = {name: 0 for name in self.names}
            print(self.names, self.name_count)
            
            return 
        
        except Exception as e:
            self.show_error_message(str(e))
            
            return
    
    
    def run_game(self):
        try:
            self.progress_bar.setValue(0)  # Reset progress bar
            self.progress_bar.setMaximum(100)
            pick_random_name(self.names, self.name_count) 
            self.progress_bar.setValue(np.random.randint(0, 100))
            return
            
        except Exception as e:
            self.show_error_message(str(e))
            self.run_button.setEnabled(True)  # Re-enable the run button if an error occurs
            return   
    
    def show_error_message(self, message):
        QMessageBox.critical(self, 'Error', message)
    
    
# Function to animate the name over the graph for a short amount of time
def animate_name(picked_name, ax):
    middle_x = len(game.names) / 2 - 0.5
    middle_y = max(game.name_count.values()) / 2

    num_stars = 15
    thetas = np.linspace(0, 2*np.pi, num_stars)
    leading_shots = max(game.name_count.values())
    
    to_loop = thetas.copy()
    # Error here collapses in some way?
    for theta in to_loop:
        
        xs = middle_x + np.cos(theta) * 0.5
        ys = middle_y + np.sin(theta) * leading_shots * 0.2

        size = random.uniform(50, 100)
        color = (random.random(), random.random(), random.random())

        ax.scatter(xs, ys, s=size, c=[color], marker='*', zorder=1000)
    
    ax.annotate(picked_name, xy=(middle_x, middle_y - 0.06 * leading_shots/2), ha='center', fontsize="16", color="Black", fontweight='bold')
    
    ax.figure.canvas.draw()
    time.sleep(1)

# Function to update the histogram
def update_histogram(ax):
    ax.clear()
    ax.bar(game.name_count.keys(), game.name_count.values())
    ax.set_ylabel('Antal Shots')
    ax.set_title('Wall of Champions')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis='x', rotation=45)
    ax.figure.canvas.draw()

# Function to pick a random name and update the count
def pick_random_name(names, name_count):
    random_name = random.choice(names)
    name_count[random_name] += 1
    update_histogram(game.figure_canvas.figure.gca())
    animate_name(random_name, game.figure_canvas.figure.gca())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    game = GameApp()
    game.show()
    sys.exit(app.exec_())
