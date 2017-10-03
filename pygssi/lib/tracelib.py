#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import platform
import matplotlib
if platform.system() == 'Darwin':
    matplotlib.use('macosx')
else:
    matplotlib.use('gtk3agg')
import numpy as np
import matplotlib.pyplot as plt


class DrawableLine:
    lock = None

    def __init__(self, ax, bg, x=None, y=None):
        self.ax = ax
        self.bg = bg
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

    def update_line(self):
        # self.ax.figure.canvas.restore_region(self.bg)
        self.line.set_data(self.dist, self.y)
        self.ax.draw_artist(self.line)
        self.ax.figure.canvas.blit(self.ax.bbox)

    def on_motion(self, event):
        if not self.press or DrawableLine.lock is not self:
            return

        self.dist.append(event.xdata)
        self.y.append(event.ydata)

        if self.line is not None:
            self.update_line()
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
            self.update_line()
            self.line.set_animated(False)

    def revert(self):
        print('revert called')
        self.dist = self.oldx[:]
        self.y = self.oldy[:]

        self.line.set_animated(True)

        self.line.set_data(self.dist, self.y)

        self.line.figure.canvas.draw()
        self.line.figure.canvas.flush_events()
        self.update_line()
        self.line.set_animated(False)

    def write(self, fn):
        with open(fn, 'w') as fout:
            for x, y in zip(self.dist, self.y):
                fout.write('{:f}, {:f}\n'.format(x, y))


class LineOnRaster(DrawableLine):

    def __init__(self):
        pass


class LineList:

    def __init__(self, ax, lat=None, lon=None, linetype=DrawableLine):
        self.ax = ax
        self.linetype = linetype
        self.lines = []
        self.connect()
        self.waiting_for_renumber = False
        self.renumber = ''
        self.bg = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)

        self.lat = lat
        self.lon = lon
        self.numtext = None

        self._new()

    def connect(self):
        self.kid = self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)

    def _new(self, event=None):
        self.lines.append(self.linetype(self.ax, bg=self.bg))
        self.lines[-1].connect()
        self.linetype.lock = self.lines[-1]
        self._update_numtext()

    def _next(self, event):
        i = self.lines.index(self.linetype.lock)
        if i == len(self.lines) - 1:
            self._new(event)
        else:
            self.linetype.lock = self.lines[i + 1]
        self._update_numtext()

    def _previous(self, event):
        i = self.lines.index(self.linetype.lock)
        if i > 0:
            self.linetype.lock = self.lines[i - 1]
        self._update_numtext()

    def _undo(self, event):
        if self.linetype.lock is not None:
            self.linetype.lock.revert()

    def _save(self, event=None):
        for i, line in enumerate(self.lines):
            fn = 'traced_line_{:d}.csv'.format(i)
            line.write(fn)

    def _update_number(self, event):
        try:
            print('Renumbering to {:d}'.format(int(self.renumber)))
            self.linetype.lock = self.lines[int(self.renumber)]
        except:
            print('Failed to pickup {:s}'.format(self.renumber))
        finally:
            self.renumber = ''
            self.waiting_for_renumber = False

    def _update_numtext(self):
        if self.numtext is None or self.linetype.lock is None:
            return
        else:
            if self.linetype.lock is not None and self.linetype.lock.line is not None:
                color = self.linetype.lock.line.get_color()
            else:
                color = 'k'
            self.numtext.set_animated(True)
            self.numtext.set_text(str(self.lines.index(self.linetype.lock)))
            self.numtext.set_color(color)
            self.ax.draw_artist(self.numtext)
            self.numtext.set_animated(False)

    def on_key(self, event):
        if event.key == 'ctrl+n':
            self._new(event)
        elif event.key == 'ctrl+r':
            self.waiting_for_renumber = True
        elif event.key == 'ctrl+x':
            self._save(event)
        elif event.key == 'enter' and self.waiting_for_renumber:
            self._update_number
        elif event.key == 'ctrl+z':
            self._undo(event)
        elif event.key.isnumeric:
            self.renumber += event.key

    def add_buttons(self, button_dict):
        if 'New' in button_dict:
            button_dict['New'].on_clicked(self._new)
        if 'Next' in button_dict:
            button_dict['Next'].on_clicked(self._next)
        if 'Previous' in button_dict:
            button_dict['Previous'].on_clicked(self._previous)
        if 'Undo' in button_dict:
            button_dict['Undo'].on_clicked(self._undo)
        if 'Save' in button_dict:
            button_dict['Save'].on_clicked(self._save)

        if 'Number' in button_dict:
            button_dict['Number'].axis('off')
            self.numtext = button_dict['Number'].text(0.1, 0.9, '0', ha='left', va='top', transform=button_dict['Number'].transAxes, fontsize=24)


class zoom_factory:

    def __init__(self, ax, base_scale=2.):
        self.in_xlim = ax.get_xlim()
        self.in_ylim = ax.get_ylim()

        self.ax = ax
        self.base_scale = base_scale

    def connect(self):
        fig = self.ax.get_figure()  # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', self.zoom_fun)

    def zoom_fun(self, event):
        # get the current x and y limits
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
        cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
        xdata = event.xdata  # get event x location
        ydata = event.ydata  # get event y location
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1. / self.base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = self.base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1.

        if (xdata - cur_xrange * scale_factor < self.in_xlim[0]) or (xdata + cur_xrange * scale_factor > self.in_xlim[1]) or (ydata - cur_yrange * scale_factor < self.in_ylim[0]) or (ydata + cur_yrange * scale_factor > self.in_ylim[1]):
            self.ax.set_xlim(self.in_xlim)
            self.ax.set_ylim(self.in_ylim)
        else:
            self.ax.set_xlim([xdata - cur_xrange * scale_factor, xdata + cur_xrange * scale_factor])
            self.ax.set_ylim([ydata - cur_yrange * scale_factor, ydata + cur_yrange * scale_factor])

        plt.draw()  # force re-draw
