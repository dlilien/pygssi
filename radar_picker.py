#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import matplotlib
matplotlib.use('tkagg')
import numpy as np
import matplotlib.pyplot as plt


class DrawableLine:
    lock = None

    def __init__(self, ax, x=None, y=None):
        self.ax = ax
        if x is not None:
            self.dist = x
            self.y = y
            self.line = self.ax.plot(self.dist, self.y)
        else:
            self.dist = []
            self.y = []
            self.line = None
        self.press = None
        self.bg = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        # if event.inaxes != self.point.axes: return
        if DrawableLine.lock is not self:
            return

        if event.xdata is None or event.ydata is None:
            return

        if event.inaxes != self.ax:
            return

        if self.line is not None:
            self.line.set_animated(True)

        self.oldx = self.dist[:]
        self.oldy = self.y[:]

        self.dist.append(event.xdata)
        self.y.append(event.ydata)

        DrawableLine.lock = self
        self.press = True

    def on_motion(self, event):
        if not self.press or DrawableLine.lock is not self:
            return

        self.dist.append(event.xdata)
        self.y.append(event.ydata)

        if self.line is not None:
            self.line.set_data(self.dist, self.y)
        else:
            self.line = self.ax.plot(self.dist, self.y)[0]

    def on_release(self, event):
        'on release we reset the press data'
        if DrawableLine.lock is not self:
            return

        self.press = None

        ind = np.argsort(self.dist)
        self.dist = [self.dist[i] for i in ind]
        self.y = [self.y[i] for i in ind]

        if self.line is None:
            self.line = self.ax.plot(self.dist, self.y)[0]
        else:
            self.line.set_data(self.dist, self.y)
            self.line.set_animated(False)
            self.line.figure.canvas.draw()

    def revert(self):
        print('revert called')
        self.dist = self.oldx[:]
        self.y = self.oldy[:]

        self.line.set_animated(True)

        self.line.set_data(self.dist, self.y)

        self.line.figure.canvas.draw()
        self.line.figure.canvas.flush_events()
        self.line.figure.canvas.draw()
        self.line.set_animated(False)


class LineList:

    def __init__(self, ax):
        self.ax = ax
        self.lines = []
        self.connect()
        self.waiting_for_renumber = False
        self.renumber = ''

    def connect(self):
        self.kid = fig.canvas.mpl_connect('key_press_event', self.on_key)

    def on_key(self, event):
        if event.key == 'ctrl+n':
            self.lines.append(DrawableLine(ax))
            self.lines[-1].connect()
            DrawableLine.lock = self.lines[-1]
        elif event.key == 'ctrl+r':
            self.waiting_for_renumber = True
        elif event.key == 'ctrl+x':
            pass
        elif event.key == 'enter':
            try:
                print('Renumbering to {:d}'.format(int(self.renumber)))
                DrawableLine.lock = self.lines[int(self.renumber)]
            except:
                print('Failed to pickup {:s}'.format(self.renumber))
            finally:
                self.renumber = ''
                self.waiting_for_renumber = False
        elif event.key == 'ctrl+z':
            print('Trying to revert')
            if DrawableLine.lock is not None:
                DrawableLine.lock.revert()
        elif event.key.isnumeric:
            self.renumber += event.key


fig = plt.figure()
plt.title('ctrl-x to save lines and exit, ctrl-n for next line\n ctrl-r, num, return, for existing number (zero-based indexing)')
ax = fig.add_subplot(111)
ax.plot(np.random.rand(10))

lines = LineList(ax)

plt.show()
