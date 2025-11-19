---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:pjrpy3] *
    language: python
    name: conda-env-pjrpy3-py
---

```python
import sys
print(sys.version)
%matplotlib inline
%run -i ~/Python/pjr3
```

```python
# emissions per decade used in the UKSM runs
em1 = np.array([0.,6.0,  18.7,  33.2,  69.9, 100.0, 250.0, 390.0])
em2 = np.array([0.,11.1,  10.6,  18.2,  42.0, 120.0, 250.0, 460.0])
em3 = np.array([0.,7.5,   9.8,  28.0,  57.0, 100.0, 280.0, 390.0])
```

```python
em1.shape
```

```python
#forcing for fixed emission estimates
emvals = np.array([0.,5.,12.5,25.,50.,100.])*4. # for 4 regions
forcvals = np.array([0.,-0.8,-1.37,-2.03,-2.95,-4.16])
```

```python
emvals
plt.plot(emvals,forcvals)
```

```python
ne = np.linspace(emvals.min(), emvals.max(), num=20, endpoint=True)
nf = np.interp(ne, emvals, forcvals)
plt.plot(ne,nf)
```

```python
f1 = np.interp(em1, emvals, forcvals)
f2 = np.interp(em2, emvals, forcvals)
f3 = np.interp(em3, emvals, forcvals)
```

```python
nyears = 70.
a1 = (f1*10).sum()/nyears
a2 = (f2*10).sum()/nyears
a3 = (f3*10).sum()/nyears
```

```python
print(a1,a2,a3)
```

```python
decade = np.array([0.,1.,2.,3.,4.,5.,6.,7.])*10+2021
plt.plot(decade,em1)
plt.plot(decade,em2)
plt.plot(decade,em3)
```

```python
plt.plot(decade,f1)
plt.plot(decade,f2)
plt.plot(decade,f3)
```

```python
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
# set labels and font size
ax.set_xlabel('Year', fontsize = 12)
ax.set_ylabel('Emissions', fontsize = 12)

#ax.plot(np.random.random(100))
#ax.plot(decade,em1,'o',label='en1')
ax.plot(decade,em1,'o',label='en1',markersize=12,fillstyle='none')
ax.plot(decade,em2,'*',label='en2')
ax.plot(decade,em3,'d',label='en3')
# change font size for x axis
#ax.xaxis.get_label().set_fontsize(20)
ax.legend()
plt.savefig('UKESM_Emissions.pdf',format='pdf',dpi=300)
plt.show()
#ax.xaxis.get_label().get_fontsize()
```

```python
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
# set labels and font size
ax.set_xlabel('Year', fontsize = 12)
ax.set_ylabel('Forcing', fontsize = 12)

#ax.plot(np.random.random(100))
ax.plot(decade,f1,'o',label='en1',markersize=12,fillstyle='none')
ax.plot(decade,f2,'*',label='en2')
ax.plot(decade,f3,'d',label='en3')

# change font size for x axis
#ax.xaxis.get_label().set_fontsize(20)
ax.legend()
#ax.xaxis.get_label().get_fontsize()
plt.savefig('UKESM_Forcing.pdf',format='pdf',dpi=300)
plt.show()
```

```python

```

```python

```
