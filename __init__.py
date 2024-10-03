
import emission

lib = emission.emission

j_fn = lib.j_fn
j_fn.__doc__ = '''Расчет эмиссии по формуле Фаулера-Нордгейма
In:  ts [K], f[V/nm], ef[eV], w[eV]
Out: J [A/m**2]
'''
j_fn.__name__ = "j_faul"

j_rd = lib.j_rds
j_rd.__doc__ = '''Расчет эмиссии по формуле Ричардсона-Дешмана-Шоттки
In:  ts [K], f[V/nm], ef[eV], w[eV]
Out: J [A/m**2]
'''
j_rd.__name__ = "j_rich"

j_c = lib.j_c
j_c.__doc__ = '''Расчет эмиссии как суммы термического и полевого потоков
In:  ts [K], f[V/nm], ef[eV], w[eV]
Out: J [A/m**2]
'''
j_c.__name__ = "j_clas"

j_q = lib.j_q
j_q.__doc__ = '''Расчет эмиссии через интегрирование решения уравнения Шредингера
In:  ts [K], f[V/nm], ef[eV], w[eV]
Out: J [A/m**2]
'''
j_q.__name__ = "j_quan"

def j_none(ts, f, ef, w):
    return 0.0

def test_emission_termal_field():
    from matplotlib.pylab import array, linspace, zeros_like
    from matplotlib.pylab import plot, show, semilogy, figure
    from matplotlib.pylab import grid, ylim, savefig, title
    tt = array([300, 500, 1000, 1500, 2000, 2500, 3000])
    ff = linspace(0.1, 10, 40)
    ef, w = 6.0, 4.5
    jj_rd = zeros_like(ff)
    jj_fn = zeros_like(ff)
    jj_c = zeros_like(ff)
    jj_q = zeros_like(ff)
    
    myplot = semilogy
    for t in tt:
        for i, f in enumerate(ff):
            jj_rd[i] = j_rd(t, f, ef, w)
            jj_fn[i] = j_fn(t, f, ef, w)
            jj_c[i]  = j_c(t, f, ef, w)
            jj_q[i]  = j_q(t, f, ef, w)

        fig = figure()
        ax = fig.add_subplot(1,1,1)
        ax.semilogy(ff, jj_rd, 'r-')
        ax.semilogy(ff, jj_fn, 'b-')
        ax.semilogy(ff, jj_c, 'g-')
        ax.semilogy(ff, jj_q, 'k-')

        ax.set_ylim([1e5, 1e14])
        ax.grid()
        ax.set_title('test_emission, t = {}'.format(t))
        ##fig.savefig('./tests/common_t_{}.png'.format(t))
        fig.show()

test_emission_termal_field()

def test_emission_holgate():
    from matplotlib.pylab import array, linspace, zeros_like
    from matplotlib.pylab import plot, show, semilogy, figure
    from matplotlib.pylab import grid, ylim, savefig, title
    
#test_emission_holgate()




##print j_fn.__doc__
##print j_rd.__doc__
##print j_c.__doc__
##print j_q.__doc__
