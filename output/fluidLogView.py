#!/usr/bin/python

import sys, pickle, pygtk, gtk, string, cairo, gobject, math
from scipy import linalg
from scipy import sparse
from gtk import *

#so, here's how we store the info in the information matrix:
# each element in it represents one timestep.
# 0: timestep #
# 1: dt
# 2: nu
# 3: kappa
# 4: initial velocity grid this timestep
# 5: initial density grid this timestep
# 6: added velocities
# 7: added densities
# 8: advected velocities
# 9: diffused velocities
# 10: q (pressures)
# 11: the gradient of q
# 12: velocity - gradient of q
# 13: q after grad q subtracted
# 14: advected densities
# 15: diffused densities

def square(x):
    return x*x
def coord2indexvel(x, y, comp, dims, z=-1):
    if len(dims) == 2:
        return comp + 2*x + 2*y*dims[0]
    else:
        return comp + 3*x + 3*y*dims[0] + 3*z*dims[1]

def coord2indexscal(x,y,dims,z=-1):
    if len(dims)==2:
        return x + y*dims[0]
    else:
        return x + y*dims[0] + z*dims[0]*dims[1]

def densterp(x,y,densities,z=0):
    if z == 0:
        grid, maxx, maxy = densities[0][0], densities[1], densities[2]
        i = int(math.floor(x-0.5))
        j = int(math.floor(y-0.5))
        s = x - i - 0.5
        r = y - j - 0.5
        pairs = [[i,j],[i+1,j],[i+1,j+1],[i,j+1]]
        
        dens_list = []
        for pair in pairs:
            tx, ty = pair[0],pair[1]
            if tx > 0 and tx < maxx and ty > 0 and ty < maxy:
                dens_list.append(grid[coord2indexscal(tx,ty,[maxx,maxy])])
            else:
                dens_list.append(0)
       # print dens_list
        #print "s:",s,"r:",r#,"x:",x,"y:",y,"i:",i,"j:",j
        q1 = dens_list[0] * (1-s) + dens_list[1] * s
        q2 = dens_list[3] * (1-s) + dens_list[2] * s
        result = q1 * (1-r) + q2 * r
        return result

class Screen(gtk.DrawingArea):

    # Draw in response to an expose-event
    __gsignals__ = { "expose-event": "override" }
    # Handle the expose-event by drawing
    def __init__(self, information):
        gtk.DrawingArea.__init__(self)
        self.set_size_request(800,800)
        self.scale_const = 25
        self.information = information
        self.timestep = 0
        self.grid_to_draw = 4 #this is the initial velocity grid
        self.dens_to_draw = 5
        self.press_to_draw = 10
        self.max_timesteps = len(information)
        self.velgrid = information[0][4]
        self.max_dens = 10.0
        self.density_xres = 50
        self.density_yres = 50
        
        
    def density_xres_handle(self,widget,box):
        self.density_xres = int(box.get_text())

    def density_yres_handle(self,widget,box):
        self.density_yres = int(box.get_text())

    def add_progress(self, progress): #we're not actually going to show() it here, just update it
        self.progress_bar = progress

    def do_expose_event(self, event):
        self.draw_grid()
        
    def timestep_entered(self, widget, box):
        tstep = int(box.get_text())
        if tstep < self.max_timesteps:
            self.timestep = tstep
            self.draw_grid()
        else:
            box.set_text("Timestep out of range.")
    def add_velstage_txtbox(self,box):
        self.velstage_box = box
    
    def add_timestep_txtbox(self,box):
        self.timestep_box = box

    def timestep_inc(self,widget,event):
        if self.timestep < (self.max_timesteps-1):
            #print "does something"
            self.timestep = self.timestep + 1
            self.timestep_box.set_text(str(self.timestep))
            self.draw_grid()
            
    def divergence_handle(self,widget,event):
        self.draw_divergence()
    
    def timestep_dec(self,widget,event):
        if self.timestep > 0:
            self.timestep = self.timestep - 1
            self.timestep_box.set_text(str(self.timestep))
            self.draw_grid()

    def dens_handle(self,widget,event):
        self.draw_density()

    def press_handle(self,widget,event):
#        print "handled
        self.press_to_draw = 10
        self.draw_pressures()

    def press_handle2(self, widget, event):
        self.press_to_draw = 13
        self.draw_pressures()

    def vel_stage_to_draw(self,widget,box):
        stage = int(box.get_text())
        acceptable = range(0,16) #[4,6,8,9,12]
        if stage in acceptable:
            if self.information[self.timestep][stage] != 'none added':
                self.grid_to_draw = stage
#                print self.grid_to_draw
                self.draw_grid()
            else:
                box.set_text("No velocity added.")
        else:
            box.set_text("Cannot draw stage " + str(stage))

    def vel_stage_inc(self,widget,event):
        current = self.grid_to_draw
        if current == 4:
            self.grid_to_draw = 8
        elif current == 8:
            self.grid_to_draw = 9
        elif current == 9:
            self.grid_to_draw = 12
        self.velstage_box.set_text(str(self.grid_to_draw))
        self.draw_grid()

    def vel_stage_dec(self,widget,event):
        current = self.grid_to_draw
        if current == 12:
            self.grid_to_draw = 9
        elif current == 9:
            self.grid_to_draw = 8
        elif current == 8:
            self.grid_to_draw = 4
        self.velstage_box.set_text(str(self.grid_to_draw))
        self.draw_grid()
    
    def foo_button(self,widget,event):
        self.draw_grid()
        
    def draw_grid(self):
        self.ctx = self.window.cairo_create()
        self.ctx.rectangle(0,0,1000,1000)
        self.ctx.set_source_rgb(1,1,1)
        self.ctx.fill()
        grid = self.information[self.timestep][self.grid_to_draw]
        idisp = self.scale_const/2 #we don't want to start drawing at the edge of the screen because that would suck
        #def coord2indexvel(x, y, comp, dims, z=-1):
        arr = grid[0]
        maxx = grid[1]
        maxy = grid[2]
        comps = grid[3]

        #let's draw a rectangle around the grid
        self.ctx.set_source_rgba(0,0,0,0.5)
        self.ctx.set_line_width(1)
        self.ctx.new_path()
        self.ctx.rectangle(idisp,idisp, maxx*self.scale_const, maxy*self.scale_const)

        if comps == 2:

            for i in range(0,maxx):
                self.ctx.move_to(idisp+i*self.scale_const,idisp)
                self.ctx.rel_line_to(0,maxy*self.scale_const)
            for j in range(0,maxy):
                self.ctx.move_to(idisp,idisp+j*self.scale_const)
                self.ctx.rel_line_to(maxx*self.scale_const,0)
            self.ctx.stroke()
            self.ctx.set_line_width(2)
            self.ctx.set_source_rgb(0,0,0)
            for i in range(0,maxx):
                for j in range(0,maxy):
                    ix = self.scale_const*i+self.scale_const
                    iy = self.scale_const*j+self.scale_const
                    fx = ix+self.scale_const/10.0*arr[coord2indexvel(i,j,0,[maxx,maxy])]
                    fy = iy+self.scale_const/10.0*arr[coord2indexvel(i,j,1,[maxx,maxy])]
                    if (fx-ix) != 0.0 or (fy-iy) != 0.0:
                        #print "attempting arrow"
                        self.draw_arrow(ix,iy,fx,fy)
#         self.ctx.move_to(0,0)
#         self.ctx.stroke()
                        


    def draw_arrow(self,ix,iy,fx,fy):
        #self.ctx = self.window.cairo_create()
        self.ctx.set_source_rgb(0,0,0)
        dv = [fx-ix,fy-iy]
        angle = math.atan(dv[0]/dv[1])
        mag = math.sqrt(square(dv[1]) + square(dv[0]))

        self.ctx.move_to(ix,iy)
        self.ctx.line_to(fx,fy)
        self.ctx.line_to(fx + (mag / 10) * (math.cos (angle + math.pi * 3 / 4)),
                     fy + (mag/10) * (math.sin (angle + math.pi * 3/4)))
        self.ctx.move_to(fx,fy)
        self.ctx.line_to(fx + (mag / 10) * (math.cos (angle - math.pi * 3 / 4)),
                     fy + (mag/10) * (math.sin (angle - math.pi * 3/4)))
        self.ctx.move_to(ix,iy)
        self.ctx.move_to(0,0)
        self.ctx.stroke()

    def draw_pressures(self):
        self.ctx = self.window.cairo_create()
#        print "in draw_pressures"
        #we'll assume 2d in here
        self.pressures = self.information[self.timestep][self.press_to_draw]
        self.grid, self.maxx, self.maxy = self.pressures[0],self.pressures[2],self.pressures[3]
        idisp = self.scale_const / 2
#        print self.grid
        for i in range(0,self.maxx):
            for j in range(0,self.maxy):
#                print i,j
                curr_press = self.grid[coord2indexscal(i,j,[self.maxx,self.maxy])]
                if curr_press > 0:
                    self.ctx.set_source_rgba(1,0,0,0.6)
                    self.ctx.arc(idisp*2+i*self.scale_const, idisp*2+j*self.scale_const, curr_press * self.scale_const / 10.0, 0, math.pi * 2)
                    self.ctx.fill()
                elif curr_press < 0:
                    self.ctx.set_source_rgba(0,0,1,0.6)
                    self.ctx.arc(i*self.scale_const+idisp*2, idisp*2+j*self.scale_const, curr_press * self.scale_const / 10.0, 0, math.pi * 2)
                    self.ctx.fill()

    def draw_density(self):
        self.dens = self.information[self.timestep][self.dens_to_draw]
        if len(self.dens) == 3:
            maxx, maxy = self.dens[1], self.dens[2]
            self.ctx = self.window.cairo_create()
            counter = 0
            idisp = self.scale_const/2.0
            maxcount = self.density_xres * self.density_yres / 100
            xrang = [self.scale_const/2, self.scale_const/2 + maxx*self.scale_const]
            yrang = [self.scale_const/2,self.scale_const/2 + maxy * self.scale_const]
            xstep = (xrang[1]-xrang[0])/self.density_xres
            ystep = (yrang[1]-yrang[0])/self.density_yres
#            print xstep, ystep
            for i in xrange(0,self.density_xres):
                for j in xrange(0,self.density_yres):
                    if counter>= maxcount:
                        counter = 0
#                        print (j+i*maxy*1.0)/(maxcount * 100)
                        #self.progress_bar.set_fraction()
                    else:
                        counter = counter + 1
                    tx = maxx *1.0 *i/ self.density_xres + maxx*1.0/self.density_xres/2.0
                    ty = maxy * 1.0*j / self.density_yres + maxy*1.0/self.density_yres/2.0
                    densi = densterp(tx,ty,self.dens)
#                    print densi
                    if densi > 0:
                        self.ctx.set_source_rgba(0.543,0.271,0.186,densi/self.max_dens)#brown
                        self.ctx.rectangle(idisp+i*xstep,idisp+j*ystep,xstep,ystep)
                        self.ctx.fill()

    def draw_divergence(self):
        print "in draw_divergence"
        self.ctx = self.window.cairo_create()
        grid = self.information[self.timestep][self.grid_to_draw]
        divmat = make_div_matrix(grid)
        divmat.tocsr()
        divergence = divmat.matvec(grid[0])
        maxx, maxy = grid[1], grid[2]
        print maxx, maxy
        idisp = self.scale_const / 2
#        print self.grid
        for i in range(0,maxx):
            for j in range(0,maxy):
                #print i,j
                curr_press = divergence[coord2indexscal(i,j,[maxx,maxy])]
                #print curr_press
                if curr_press > 0:
                    self.ctx.set_source_rgba(1,0,0,0.6)
                    self.ctx.arc(idisp*2+i*self.scale_const, idisp*2+j*self.scale_const, curr_press * self.scale_const / 10.0, 0, math.pi * 2)
                    self.ctx.fill()
                elif curr_press < 0:
                    self.ctx.set_source_rgba(0,0,1,0.6)
                    self.ctx.arc(i*self.scale_const+idisp*2, idisp*2+j*self.scale_const, abs(curr_press) * self.scale_const / 10.0, 0, math.pi * 2)
                    self.ctx.fill()
        

def pdbwrap(f):
    """
    Wrap a closure to cause it to dump into the debugger when there's an exception.
    I use this in GTK callbacks, for instance:
        self.connect('key_press_event', pdbwrap(self.key_press_cb))
    Otherwise, the pygtk code just prints the error and continues, making it hard to debug
    """
    def fdebug(*a, **kw):
        try:
            return f(*a, **kw)
        except Exception:
            type, value, tb = sys.exc_info()
            traceback.print_exc(file=sys.stderr)
            if sys.stdin.isatty():
                pdb.post_mortem(tb)
            else:
                sys.exit(1)
    return fdebug

def make_div_matrix(grid):
    if grid[-1]==2:
        maxx, maxy = grid[1],grid[2]
        divmat = sparse.lil_matrix((maxx*maxy,2* maxx*maxy))
        for i in range(0,maxx):
            for j in range(0,maxy):
                row = i + j*maxx
                x1,x2 = [i-1,j], [i,j]
                if x1[0] > 0 and x1[0] < maxx and x1[1]>0 and x1[1] < maxy:
                    divmat[row,coord2indexvel(x1[0],x1[1],0,[maxx,maxy])] = -1
                if x2[0] > 0 and x2[0] < maxx and x2[1]>0 and x2[1] < maxy:
                    divmat[row,coord2indexvel(x2[0],x2[1],0,[maxx,maxy])] = 1
                y1,y2 = [i,j-1],[i,j]
                if y1[0]>0 and y1[0] < maxx and y1[1]>0 and y1[1] < maxy:
                    divmat[row,coord2indexvel(y1[0],y1[1],1,[maxx,maxy])] = -1
                if y2[0]>0 and y2[0] < maxx and y2[1]>0 and y2[1] < maxy:
                    divmat[row,coord2indexvel(y2[0],y2[1],1,[maxx,maxy])] = 1
        return divmat

    
def main():
    
    fname = sys.argv[1]
    fstream = open(fname, 'r')
    information = pickle.load(fstream)
    fstream.close()
    
#    print information[0][6]
    
    # for elt in grid[0]:
#         print elt
    
    window = gtk.Window()
    window.connect("delete-event", gtk.main_quit)
    widget = Screen(information)
    
    bigbox = gtk.HBox()
    
    buttonbox = gtk.VBox()
    button = gtk.Button("Redraw Grid")
    handler = button.connect('button_press_event', widget.foo_button)

    

    textlabel = gtk.Label("\nTimestep")
    timestep = gtk.Entry(40)
    timestep.set_text("0")
    timestep_entry = timestep.connect('activate', widget.timestep_entered, timestep)

    timestepButtonBox = gtk.HBox()
    timestepPlusButton = gtk.Button("+")
    timestepIncHandler = timestepPlusButton.connect('button_press_event',widget.timestep_inc)
    timestepMinusButton = gtk.Button("-")
    timestepDecHandler = timestepMinusButton.connect('button_press_event', widget.timestep_dec)
    timestepButtonBox.pack_start(timestepMinusButton)
    timestepButtonBox.pack_start(timestepPlusButton)

    widget.add_timestep_txtbox(timestep)
    stagelabel = gtk.Label("\nVelocity Stage")
    stage = gtk.Entry(40)
    stage.set_text("4")
    stage_entry = stage.connect('activate',widget.vel_stage_to_draw, stage)
    widget.add_velstage_txtbox(stage)
    stage_inc = gtk.Button("+")
    stage_inc_handle = stage_inc.connect('button_press_event', widget.vel_stage_inc)
    stage_dec = gtk.Button("-")
    stage_dec_handle = stage_dec.connect('button_press_event', widget.vel_stage_dec)
    stage_box = gtk.HBox()
    stage_box.pack_start(stage_dec)
    stage_box.pack_start(stage_inc)
    newline = gtk.Label(" ")

    density_butt = gtk.Button("Draw density")
    density_butt_handle = density_butt.connect('button_press_event', widget.dens_handle)
    
    density_res_box = gtk.HBox()
    density_xres = gtk.Entry(4)
    density_yres = gtk.Entry(4)
    density_res_box.pack_start(density_xres,False)
    density_res_box.pack_start(density_yres,False)
    density_xres_handler = density_xres.connect('activate',widget.density_xres_handle,density_xres)
    density_yres_handler = density_yres.connect('activate',widget.density_yres_handle,density_yres)

    density_progress = gtk.ProgressBar()
    widget.add_progress(density_progress)

    presslabel = gtk.Label("\nPressure")
    pressbox = gtk.HBox()
    pressure_butt = gtk.Button("Before")
    press_handle = pressure_butt.connect('button_press_event', widget.press_handle)
    pressure_butt2 = gtk.Button("After")
    press_handle = pressure_butt.connect('button_press_event', widget.press_handle2)
    pressbox.pack_start(pressure_butt)
    pressbox.pack_start(pressure_butt2)

    
    divergence = gtk.Button("Divergence")
    divergence_handle = divergence.connect('button_press_event', widget.divergence_handle)


    
    buttonbox.pack_start(button,False)
    buttonbox.pack_start(textlabel,False)
    buttonbox.pack_start(timestep,False)
    buttonbox.pack_start(timestepButtonBox, False)
    buttonbox.pack_start(stagelabel,False)
    buttonbox.pack_start(stage,False)
    buttonbox.pack_start(stage_box, False)
    buttonbox.pack_start(newline, False)
    buttonbox.pack_start(density_butt,False, padding=3)
    buttonbox.pack_start(density_res_box, False)
    buttonbox.pack_start(density_progress,False)
    buttonbox.pack_start(presslabel,False)
    buttonbox.pack_start(pressbox, False)
    buttonbox.pack_start(divergence,False)
#    buttonbox.pack_start(pressure_butt, False)
    bigbox.pack_end(buttonbox, False, True, 5)
    bigbox.pack_end(widget, True, True, 0)
    buttonbox.show()
    stage_inc.show()
    stage_box.show()
    density_res_box.show()
    density_xres.show()
    density_yres.show()
    newline.show()
    stage_dec.show()
    divergence.show()
    density_butt.show()
    pressure_butt2.show()
    pressbox.show()
    density_progress.show()
    timestepPlusButton.show()
    timestepMinusButton.show()
    timestepButtonBox.show()
    pressure_butt.show()
    presslabel.show()
    stage.show()
    stagelabel.show()
    timestep.show()
    button.show()
    widget.show()
    textlabel.show()
    window.add(bigbox)
    bigbox.show()
    window.present()
    gtk.main()

if __name__ == "__main__":
    #print grid
    
#    pdbwrap(main)()
    main()
