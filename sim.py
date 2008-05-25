#!/usr/bin/python

from numpy import allclose, arange, eye, linalg, random, ones, zeros
from scipy import linalg, sparse

import cairo, gtk, math

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys

def square(x):
    return x*x

scale_const = 25

max_dens = 10.0

grid = 0

# Some api in the chain is translating the keystrokes to this octal string
# so instead of saying: ESCAPE = 27, we use the following.
ESCAPE = '\033'
SPACE = '\040'

# Number of the glut window.
window = 0

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


def append_meta(grid, x, y, comps=1, zee=0):
    if zee>0:
        return [grid,x,y,zee,comps]
    else:
        return [grid,x,y,comps]

def build_force_grid(xrange,yrange,x,y,maxx,maxy):
    arr = zeros(2*maxx*maxy)
    for i in xrange:
        for j in yrange:
            xc = coord2indexvel(i,j,0,[maxx,maxy])
            yc = coord2indexvel(i,j,1,[maxx,maxy])
            arr[xc] = x
            arr[yc] = y
    return arr

def build_gravity(x,y,z=0):
    xgrav = 0.1
    ygrav = -5.0
    return build_force_grid(range(0,x),range(0,y),xgrav,ygrav,x,y)

def build_densities(x,y,z=0):
    arr = zeros(x*y)
    for i in range((x-1)/3+1,2*(x-1)/3+1):
        for j in range(2*(x-1)/3+1,y):
            arr[coord2indexscal(i,j,(x,y))] = max_dens + 5
    return [[arr,5],x,y]

def sum_density(densities):
    dens = 0
    grid = densities[0][0]
    if len(densities) == 3:
        for i in xrange(0,densities[1]):
            for j in xrange(0,densities[2]):
                dens += grid[coord2indexscal(i,j,[densities[1],densities[2]])]
    return dens
                
def build_vel(x,y,z=0):
    if z > 0:
        arr = zeros(3*x*y*z)
    else:
        arr = zeros(2*x*y)
    if z==0:
        for i in range((x-1)/3+1,2*(x-1)/3+1):
            for j in range(2*(y-1)/3+1,y):
                xc = coord2indexvel(i,j,0,[x,y])
                yc = coord2indexvel(i,j,1,[x,y])
                arr[xc] = 0.0
                arr[yc] = -10.0
    else:
        for i in range((x-1)/3+1,2*(x-1)/3+1):
            for j in range(0,(y-1)/3+1):
                for k in range(0,(z-1)/3+1):
                    xc = coord2indexvel(i,j,0,[x,y,z],k)
                    yc = coord2indexvel(i,j,1,[x,y,z],k)
                    zc = coord2indexvel(i,j,2,[x,y,z],k)
                    arr[xc] = 5.0
                    arr[yc] = 1.0
                    arr[zc] = 0.5
    if z>0:
        return append_meta(arr,x,y,zee=z,comps=3)
    else:
        return append_meta(arr,x,y,2)

def advect_grid(grid,dt):
    if grid[-1] == 2: 
        arr = grid[0]
        maxx = grid[1]
        maxy = grid[2]
        newgrid = zeros(2*maxx*maxy)
        for i in range(0,maxx):
            for j in range(0,maxy):
                xvel = arr[coord2indexvel(i,j,0,[maxx,maxy])]
                yvel = arr[coord2indexvel(i,j,1,[maxx,maxy])]
                newx = i - xvel * dt
                newy = j - yvel * dt

                vel = linterp(newx, newy, grid)
                if vel == 0:
                    vel = [0.0,0.0]
                newgrid[coord2indexvel(i,j,0,[maxx,maxy])] = vel[0]
                newgrid[coord2indexvel(i,j,1,[maxx,maxy])] = vel[1]
        return [newgrid,maxx,maxy,2]

#def coord2indexvel(x, y, comp, dims, z=-1):
def linterp(x,y,grid,z=0):
    if z==0:
        arr = grid[0]
        maxx = grid[1]+1
        maxy = grid[2]+1
        arrmaxx = grid[1]
        arrmaxy = grid[2]
        retval = []
        if x<0 or x>maxx or y<0 or y>maxy:
            return 0
        else:
            corner = [int(math.floor(x)),int(math.floor(y-0.5))]
            a,b = corner[0],corner[1]
            r=y-corner[1]
            s=x-corner[0]
            R = r-0.5
            S = s

            pairs = [[a-1,b],[a,b],[a,b+1],[a-1,b+1]]
            velocities = []
            for pair in pairs:
                xt,yt = pair[0], pair[1]
                if xt < 0 or yt < 0 or xt >= arrmaxx or yt >= arrmaxy:
                    velocities.append(0)
                else:
                    velocities.append(arr[coord2indexvel(xt,yt,0,[arrmaxx,arrmaxy])])
            E = velocities[0] * (1-S) + velocities[1] * S
            D = velocities[3] * (1-S) + velocities[2] * S
            F = R * D + (1-R) * E
            retval.append(F)
            A = int(math.floor(x-0.5))
            B = int(math.floor(y))
            del velocities[:]
            R = x-A-0.5
            S = y-B
            pairs = [[A,B-1],[A+1,B-1],[A+1,B],[A,B]]
            for pair in pairs:
                xt = pair[0]
                yt = pair[1]
                if xt < 0 or yt < 0 or xt >= arrmaxx or yt >= arrmaxy:
                    velocities.append(0)
                else:
                    velocities.append(arr[coord2indexvel(xt,yt,1,[arrmaxx,arrmaxy])])
            E = (1-R) * velocities[0] + R * velocities[1]
            D = (1-R) * velocities[3] + R * velocities[2]
            F = (1-S) * E + S * D
            retval.append(F)
            return retval


def make_vector_lapl_mat(grid):
    if grid[-1] == 2:
        maxx, maxy = grid[1], grid[2]
        vect = sparse.lil_matrix((2*maxx*maxy, 2*maxx*maxy))
        for i in range(0,maxx):
            for j in range(0,maxy):
                row = 2 * i + 2 * j * maxx
                pairs = [[i-1,j],[i,j-1],[i,j],[i,j+1],[i+1,j]]
                for pair in pairs:
                    tmpx, tmpy = pair[0], pair[1]
                    if tmpx > 0 and tmpx < maxx and tmpy > 0 and tmpy < maxy:
                        if tmpx == i and tmpy == j:
                            vect[row, coord2indexvel(tmpx, tmpy, 0, [maxx, maxy])] = -4
                            vect[row+1, coord2indexvel(tmpx, tmpy, 1, [maxx, maxy])] = -4
                        else:
                            vect[row,coord2indexvel(tmpx, tmpy, 0, [maxx, maxy])] = 1
                            vect[row+1,coord2indexvel(tmpx, tmpy, 1, [maxx, maxy])] = 1
        return vect

def make_diffuse_matrix(grid, nu, dt):
    laplacian = make_vector_lapl_mat(grid)
    maxx, maxy = grid[1], grid[2]
    tmp_mat = sparse.lil_matrix((2*maxx*maxy, 2*maxx*maxy))
    tmp_mat.setdiag(ones(2*maxx*maxy))
    laplacian *= nu * dt
    tmp_mat -= laplacian
    return tmp_mat

def make_dens_diff_matrix(grid,nu,dt):
    laplacian = make_scalar_lapl(grid)
    maxx, maxy = grid[1], grid[2]
    tmp_mat = sparse.lil_matrix((maxx*maxy, maxx*maxy))
    tmp_mat.setdiag(ones(maxx*maxy))
    laplacian *= nu * dt
    tmp_mat -= laplacian
    return tmp_mat

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

def make_scalar_lapl(grid):
    if grid[-1]==2:
        maxx, maxy = grid[1], grid[2]
        scalapl = sparse.lil_matrix((maxx*maxy, maxx*maxy))
        for i in range(0,maxx):
            for j in range(0,maxy):
                row = i+j*maxx
                pairs = [[i-1,j],[i+1,j],[i,j-1],[i,j+1],[i,j]]
                for pair in pairs:
                    tmpx, tmpy = pair[0],pair[1]
                    if tmpx >= 0 and tmpx < maxx and tmpy >= 0 and tmpy < maxy:
                        if tmpx == i and tmpy == j:
                            scalapl[row, coord2indexscal(tmpx, tmpy,[maxx, maxy])] = -4
                        else:
                            scalapl[row,coord2indexscal(tmpx, tmpy,[maxx, maxy])] = 1
        return scalapl

def make_grad_mat(grid):
    if grid[-1]==2:
        maxx, maxy = grid[1], grid[2]
        grad = sparse.lil_matrix((2*maxy*maxx, maxx*maxy))
        for i in range(0,maxx):
            for j in range(0,maxy):
                xrow = 2*i + 2*j*maxx
                if i >=0 and j >= 0 and i < maxx and j < maxy:
                    grad[xrow, coord2indexscal(i,j,[maxx,maxy])] = -1
                if (i+1) >=0 and j >= 0 and (i+1) < maxx and j < maxy:
                    grad[xrow,coord2indexscal(i+1,j,[maxx,maxy])] = 1
                yrow = xrow+1
                if i >=0 and j >= 0 and i < maxx and j < maxy:
                    grad[yrow,coord2indexscal(i,j,[maxx,maxy])] = -1
                if i >=0 and (j+1) >= 0 and i < maxx and (j+1) < maxy:
                    grad[yrow, coord2indexscal(i,j+1,[maxx,maxy])] = 1
        return grad


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

def advect_densities(densities,velocity,dt):
    if velocity[-1] == 2: 
        velarr = velocity[0]
        maxx = velocity[1]
        maxy = velocity[2]
        newgrid = zeros(maxx*maxy)
        for i in xrange(0,maxx):
            for j in xrange(0,maxy):
                xvel = velarr[coord2indexvel(i,j,0,[maxx,maxy])]
                yvel = velarr[coord2indexvel(i,j,1,[maxx,maxy])]
                newx = i - xvel * dt
                newy = j - yvel * dt
                newgrid[coord2indexscal(i,j,[maxx,maxy])] = densterp(newx,newy,densities)
        return [[newgrid,5],maxx,maxy]

#
#
##
## DRAWING CRAP FOLLOWING              ******************************************************************************
## THIS LINE                           ********************************************************************************
##                                     ********************************************************************************
#
#
def draw_arrow(ix,iy,fx,fy,surf):
    dv = [fx-ix,fy-iy]
    angle = math.atan(dv[0]/dv[1])
    mag = math.sqrt(square(dv[1]) + square(dv[0]))
    
    surf.move_to(ix,iy)
    surf.line_to(fx,fy)
    surf.line_to(fx + (mag / 10) * (math.cos (angle + math.pi * 3 / 4)),
                 fy + (mag/10) * (math.sin (angle + math.pi * 3/4)))
    surf.move_to(fx,fy)
    surf.line_to(fx + (mag / 10) * (math.cos (angle - math.pi * 3 / 4)),
                 fy + (mag/10) * (math.sin (angle - math.pi * 3/4)))
    surf.move_to(ix,iy)


#so we're assuming a 2d grid. since we're not using cairo for 3d grids, 
#this is perfectly reasonable
def draw_grid(surface,grid):
    ctx = surface
    ctx.rectangle(0,0,1000,1000)
    ctx.set_source_rgb(1,1,1)
    ctx.fill()
    
    idisp = scale_const/2 #we don't want to start drawing at the edge of the screen because that would suck
    #def coord2indexvel(x, y, comp, dims, z=-1):
    arr = grid[0]
    maxx = grid[1]
    maxy = grid[2]
    comps = grid[3]

    #let's draw a rectangle around the grid
    ctx.set_source_rgba(0,0,0,0.5)
    ctx.set_line_width(1)
    ctx.new_path()
    ctx.rectangle(idisp,idisp, maxx*scale_const,maxy*scale_const)

    if comps == 2:

        for i in range(0,maxx):
            ctx.move_to(idisp+i*scale_const,idisp)
            ctx.rel_line_to(0,maxy*scale_const)
        for j in range(0,maxy):
            ctx.move_to(idisp,idisp+j*scale_const)
            ctx.rel_line_to(maxx*scale_const,0)
        ctx.stroke()
        ctx.set_line_width(2)
        ctx.set_source_rgb(0,0,0)
        for i in range(0,maxx):
            for j in range(0,maxy):
                ix = scale_const*i+scale_const
                iy = scale_const*j+scale_const
                fx = ix+scale_const/10*arr[coord2indexvel(i,j,0,[maxx,maxy])]
                fy = iy+scale_const/10*arr[coord2indexvel(i,j,1,[maxx,maxy])]
                if (fx-ix) != 0.0 and (fy-iy) != 0.0:
                    draw_arrow(ix,iy,fx,fy,ctx)
    ctx.move_to(0,0)
    ctx.stroke()

def draw_pressures(surface, pressures):
    ctx = surface
    #we'll assume 2d in here
    grid, maxx, maxy = pressures[0][0],pressures[1],pressures[2]
    idisp = scale_const / 2
    for i in range(0,maxx):
        for j in range(0,maxy):
            curr_press = grid[coord2indexscal(i,j,[maxx,maxy])]
            if curr_press > 0:
                ctx.set_source_rgba(1,0,0,0.6)
                ctx.arc(idisp*2+i*scale_const, idisp*2+j*scale_const, curr_press * scale_const / 10.0, 0, math.pi * 2)
                ctx.fill()
            elif curr_press < 0:
                ctx.set_source_rgba(0,0,1,0.6)
                ctx.arc(i*scale_const+idisp*2, idisp*2+j*scale_const, curr_press * scale_const / 10.0, 0, math.pi * 2)
                ctx.fill()


def button_press_event(widget,event,grid):
#    print 'pressed, location:', event.x, event.y
    tx = event.x
    ty = event.y
    idisp = scale_const/2
    rx = (tx-scale_const/2)/scale_const
    ry = (ty-scale_const/2)/scale_const
 #   print rx,ry
#    print densterp(rx,ry,grid)
    return True


def draw_density(surface,density):
    if len(density) == 3:
        maxx, maxy = density[1], density[2]
        ctx = surface
        counter = 0
        maxcount = maxx*maxy*scale_const*scale_const / 100
        for ipix in xrange(scale_const/2,scale_const/2 + maxx * scale_const):
            for jpix in xrange(scale_const/2,scale_const/2 + maxy * scale_const):
                if counter >= maxcount:
                    counter = 0
                    print (jpix*1.0+ipix*maxy*scale_const)/(maxx*maxy*scale_const*scale_const)*100,"%" 
                else:                    
                    counter = counter + 1
                tx = (ipix-scale_const/2.0)/scale_const
                ty = (jpix-scale_const/2.0)/scale_const
                dens = densterp(tx,ty,density)
                ctx.set_source_rgba(0.543,0.271,0.186,dens/max_dens)#brown
                ctx.rectangle(ipix,jpix,1,1)
                ctx.fill()


##
##
##
###
#### BOILERPLATE
####
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##
##
##

def InitGL(Width, Height):
    glClearColor(0.0, 0.0, 0.0, 0.0)
    glClearDepth(1.0)
    glDepthFunc(GL_LESS)
    glEnable(GL_DEPTH_TEST)
    glShadeModel(GL_SMOOTH)
    glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
    glEnable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)

    glMatrixMode(GL_MODELVIEW)

def ReSizeGLScene(width, height):
    if height == 0:
        height = 1

    glViewport(0,0,width,height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, float(width)/float(height), 0.1, 100.0)
    
    glMatrixMode(GL_MODELVIEW)

# The main drawing function. 
def DrawGLScene():
    # Clear The Screen And The Depth Buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()					# Reset The View 

    # Move Left 1.5 units and into the screen 6.0 units.
    glTranslatef(0.0, 0.0, -8.0)

    glColor3f(1.0,0.0,0.0)
    #glRectf(-3.0,2.0,3.0,-2.0)
    if grid[-1] == 2:
        glBegin(GL_LINE_STRIP)
        glVertex3f(-3.0,3.0,0.0)
        glVertex3f(3.0,3.0,0.0)
        glVertex3f(3.0,-3.0,0.0)
        glVertex3f(-3.0,-3.0,0.0)
        glVertex3f(-3.0,3.0,0.0)
        glEnd()
        
        arr, maxx, maxy = grid[0],grid[1],grid[2]
        #print maxx, maxy
        xdisp = 6.0 / (maxx+1)
        ydisp = 6.0 / (maxy+1)
        
        #print xdisp, ydisp
        
        glTranslatef(-3.0+xdisp,3.0-ydisp,0.0)
        for j in xrange(0,maxy):
            for i in xrange(0,maxx):
                xv = arr[coord2indexvel(i,j,0,[maxx,maxy])] / scale_const / 3
                yv = -1.0 * arr[coord2indexvel(i,j,1,[maxx,maxy])] / scale_const / 3
                #print xv,yv
                glBegin(GL_LINES)
                glVertex3f(0.0,0.0,0.0)
                glVertex3f(xv,yv,0.0)
                glEnd()
                
                #now we move the initial point (make sure to take care of y-edge-case (need a new line)
                if i == (maxx-1):
                    #print i
                    glTranslatef((-30.0)*xdisp+xdisp, -1 * ydisp,0.0)
                else:
                    glTranslatef(xdisp, 0.0,0.0)
        glTranslatef(-xdisp-0.05,ydisp*31+0.1,0.0)
        
        if draw_density == 1:
            dx = 6.0 / dens_squares_x
            dy = 6.0 / dens_squares_y
            counter = 0
            for j in xrange(0,dens_squares_y):
                for i in xrange(0,dens_squares_x):
                    if counter == 500:
                        print (i+j*dens_squares_y*1.0)/(dens_squares_x*dens_squares_y)*100,"%"
                        counter = 0
                    else:
                        counter += 1
                    tx = 1.0 * i / dens_squares_x * maxx
                    ty = 1.0 * j / dens_squares_y * maxy
                    curr_dens = densterp(tx,ty,density)
                    

                    glColor4f(0.543,0.271,0.186,(1.0 * curr_dens)/max_dens)
                    #glBlendFunc(GL_SRC_ALPHA, GL_ONE)
                    glBegin(GL_QUADS)
                    glVertex3f(-dx,dy,0.0)
                    glVertex3f(dx,dy,0.0)
                    glVertex3f(dx,-dy,0.0)
                    glVertex3f(-dx,-dy,0.0)
                    glEnd()
#                     if curr_dens != 0.0:
#                         glColor3f(0.0,1.0,0.0)
#                         glBegin(GL_LINES)
#                         glVertex3f(0.0,0.0,0.0)
#                         glVertex3f(1.0,1.0,0.0)
#                         glEnd()
                    
                    if i == (dens_squares_x - 1):
                        glTranslatef(-6.0+dx,-dy,0.0)
                    else:
                        glTranslatef(dx,0.0,0.0)
                    
    glutSwapBuffers()

# The function called whenever a key is pressed. Note the use of Python tuples to pass in: (key, x, y)  
def keyPressed(*args):
	# If escape is pressed, kill everything.
    if args[0] == ESCAPE:
        sys.exit()
    if args[0] == SPACE:
        update_grids()
    if args[0] == '\141': #a
        add_source()
    if args[0] == '\163':
        density[0][0] += source_dens[0][0]
def main():
    global window
    # For now we just pass glutInit one empty argument. I wasn't sure what should or could be passed in (tuple, list, ...)
    # Once I find out the right stuff based on reading the PyOpenGL source, I'll address this.
    glutInit(sys.argv)

    # Select type of Display mode:   
    #  Double buffer 
    #  RGBA color
    # Alpha components supported 
    # Depth buffer
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    
    # get a 640 x 480 window 
    glutInitWindowSize(800, 800)
    
    # the window starts at the upper left corner of the screen 
    glutInitWindowPosition(0, 0)
    
    # Okay, like the C version we retain the window id to use when closing, but for those of you new
    # to Python (like myself), remember this assignment would make the variable local and not global
    # if it weren't for the global declaration at the start of main.
    window = glutCreateWindow("Smoke")

    # Register the drawing function with glut, BUT in Python land, at least using PyOpenGL, we need to
    # set the function pointer and invoke a function to actually register the callback, otherwise it
    # would be very much like the C version of the code.	
    glutDisplayFunc(DrawGLScene)
    
    # Uncomment this line to get full screen.
    #glutFullScreen()
    # When we ar doing nothing, redraw the scene.
    #glutIdleFunc(DrawGLScene)
    
    # Register the function called when our window is resized.
    glutReshapeFunc(ReSizeGLScene)
    
    # Register the function called when the keyboard is pressed.  
    glutKeyboardFunc(keyPressed)
    # Initialize our window 
    InitGL(800,800)
    # Start Evet Processing Engine	
    glutMainLoop()





 #in pixels
dens_squares_x = 90
dens_squares_y = 90

draw_density = 1
grid = build_vel(30,30)

source = build_vel(30,30)
source_dens = build_densities(30,30)
max_x = grid[1]
max_y = grid[2] 
dens_diff = make_dens_diff_matrix(grid,0.5,0.1)
dens_diff = dens_diff.tocsr()
density = build_densities(max_x,max_y)

dens_before = sum_density(density)
grav = build_gravity(max_x,max_y)
diff_mat = make_diffuse_matrix(grid,0.5,0.1)
diff_mat = diff_mat.tocsr()
div_mat = make_div_matrix(grid)
div_mat = div_mat.tocsr()
scalp = make_scalar_lapl(grid)
scalp = scalp.tocsr()
grad_mat = make_grad_mat(grid)
grad_mat = grad_mat.tocsr()

def add_source():
    global grid
    grid[0] += source[0]

def update_grids():
    global grid, density
    grid = advect_grid(grid,0.1)
    print "Advected grid." 
    grid[0] = linalg.cg(diff_mat, grid[0])[0]
    print "Finished diffusion."
    pressures = [linalg.cg(scalp,div_mat.matvec(grid[0])), max_x, max_y]
    print "Finished calculating pressure."
    velsub = grad_mat.matvec(pressures[0][0])
    print "Finished calculating pressure gradient"
    grid[0] -= velsub
    print "Finished subtracting pressure gradient."
    pressures = [linalg.cg(scalp,div_mat.matvec(grid[0])), max_x, max_y]
    print "Finished calculating pressures for new velocity field."
    density[0][0]=linalg.cg(dens_diff,density[0][0])[0]
    print "Finished diffusing density."
    density = advect_densities(density,grid,0.1)
    print "Finished advecting density."

if __name__ == '__main__':
    # Print message to console, and kick off the main to get it rolling.
    print "Hit ESC key to quit."
    #global grid
    
    # for z in range(0,10):
#         print "step:",z
#         # grid[0] -= grav*0.05
#         #        print "Added gravity."
        
#         # if z < 4:
#         # grid[0] = grid[0] + source[0]
#         #             print "Added source force"
#         update_grids()
        
    print "Drawing density."
    main()



    
