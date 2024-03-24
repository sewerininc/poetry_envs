import sys
import random
import time
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel
from PyQt5.QtCore import Qt, QTimer
from IPython import display

# Function to animate the name over the graph for a short amount of time
def animate_name(picked_name):
    # Get the middle position of the graph
    middle_x = len(names) / 2
    middle_y = max(name_count.values()) / 2

    # Generate random number of stars
    num_stars = random.randint(5, 10)
    
    for _ in range(num_stars):
        # Generate random position around the middle
        x = middle_x + random.uniform(-0.2, 0.2)
        y = middle_y + random.uniform(-0.2, 0.2)

        # Generate random size and color
        size = random.uniform(50, 100)
        color = (random.random(), random.random(), random.random())

        ax.scatter(x, y, s=size, c=[color], marker='*')

    display.clear_output(wait=True)
    display.display(plt.gcf())
    time.sleep(1)  # Adjust the duration as needed

# Function to update the histogram
def update_histogram():
    ax.clear()
    ax.bar(name_count.keys(), name_count.values())
    ax.set_xlabel('Names')
    ax.set_ylabel('Counts')
    ax.set_title('Name Counts')
    ax.tick_params(axis='x', rotation=45)
    display.clear_output(wait=True)
    display.display(plt.gcf())

# Function to pick a random name and update the count
def pick_random_name():
    # Pick a random name
    random_name = random.choice(names)
    # Update the count for the picked name
    name_count[random_name] += 1
    # Update the histogram
    update_histogram()
    # Animate the name
    animate_name(random_name)

class GameApp(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle('Name Picker Game')
        self.setGeometry(100, 100, 400, 200)

        self.name_label = QLabel('Enter a name:')
        self.name_input = QLineEdit()
        self.run_button = QPushButton('Run')

        self.run_button.clicked.connect(self.run_game)

        vbox = QVBoxLayout()
        vbox.addWidget(self.name_label)
        vbox.addWidget(self.name_input)
        vbox.addWidget(self.run_button)

        self.setLayout(vbox)

    def run_game(self):
        name = self.name_input.text()
        if name.strip():
            names.append(name.strip())
            pick_random_name()
            self.name_input.clear()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    game = GameApp()
    game.show()
    sys.exit(app.exec_())

