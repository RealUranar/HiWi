import eel

# name of folder where the html, css, js, image files are located
eel.init('templates')

@eel.expose
def demo(x):
    return x**2

# 1000 is width of window and 600 is the height
eel.start('index.html', size=(1000, 600))

