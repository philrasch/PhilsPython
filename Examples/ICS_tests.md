---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.8
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
# imports
from icalendar import Calendar, Event, vCalAddress, vText
from datetime import datetime
from pathlib import Path
import os
import pytz
 
e = open('/home/phil/Desktop/WindowsDesktop/UWevent.ics', 'rb')
ical = Calendar.from_ical(e.read())
ocal = ical
ocal['attendee'] = ['MAILTO:maxm@mxm.dk','MAILTO:test@example.com']

for k,v in ical.items():
    print('YYY',k,v) 
for component in ocal.walk():
    print('xxx',component)
    print('component.name is ',component.name)
    print('keys',component.keys())jjfjfjfj]]
    aw
    for k,v in component.items():
        print('YY2',k,v)
    #for subcomp in component.walk():
    #    print('aaa',subcomp)
    #help(component)
    if component.name == "VEVENT":
        #component['attendee'] = ['MAILTO:maxm@mxm.dk','MAILTO:test@example.com']
        print("name",component.get("name"))
        print("desc",component.get("description"))
        print("org",component.get("organizer"))
        print("loc",component.get("location"))
        print("dtstart",component.decoded("dtstart"))
        print("dtend",component.decoded("dtend"))
    print('zzz', component)
e.close()

```

```python

```
