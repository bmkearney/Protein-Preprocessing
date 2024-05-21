#
# Add command "closest" that finds the closest pair of atoms given two
# sets of atoms.
#
# open 7LUP
# closest /A to /B
#
#  Minimum distance 2.27 between /A LYS 45 NZ and /B ASP 522 OD1
#

def humanize(session, show = False):
    msg = 'Called the script!'
    session.logger.status(msg, log = True)
    h_seq="MDTFSTKSLALQAQKKLLSKMASKAVVAVLVDDTSSEVLDELYRATREFTRSRKEAQKMLKNLVKVALKLGLLLRGDQLGGEELALLRRFRHRARCLAMTAVSFHQVDFTFDRRVLAAGLLECRDLLHQAVGPHLTAKSHGRINHVFGHLADCDFLAALYGPAEPYRSHLRRICEGLGRMLDEGSL"
    o_seq="MDSFSTKNLALQAQKKLMSKMATKTVANLFIDDTSSEVLDELYRVTKEYTRNRKEAQKIIKNLIKMVVKLGVLYRNGQFNNEELALVERFRKKVHTLAMTAVSFYQIDFTFDRRVMSNLLNDCRELLHQAINRHLTAKSHARINHVFNHFADCDFLATLYGPSEVYRGHLQKICEGVNKMLDEGNL"
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    aanum=1
    da={value: key for key, value in d.items()}
    from chimerax.core.commands import run
    for i in h_seq:
        code=da[i]
        if (h_seq[aanum-1]==o_seq[aanum-1]):
            msg = "position %d aa is same"%aanum
            session.logger.status(msg, log=True)
        else:     
            try:
                msg = "position: %d Human is "%aanum+h_seq[aanum-1]+ " Fish is "+o_seq[aanum-1]
                session.logger.status(msg, log = True)
                run(session, 'swapaa :'+str(aanum)+' '+code)
            except:
                msg = "aa doesn't exist"
                session.logger.status(msg, log = True)
        aanum+=1



    

def report_closest(session, d, a1, a2, show, max_dist):
    if d is None:
        msg = 'No atom pairs within distance %.3g of each other' % max_dist
    else:
        msg = 'Minimum distance %.2f between %s and %s' % (d, str(a1), str(a2))
    session.logger.status(msg, log = True)

    if show and a1 and a2:
        from chimerax.core.commands import run
        run(session, 'distance %s %s' % (a1.atomspec, a2.atomspec))
        run(session, 'view %s %s' % (a1.atomspec, a2.atomspec))

def closest_slow(session, atoms, to_atoms, max_dist = 10, show = False):
    '''Loop through every pair of atoms in Python.  This is slow.'''
    dmin = amin1 = amin2 = None
    from chimerax.geometry import distance
    for a1 in atoms:
        for a2 in to_atoms:
            d = distance(a1.scene_coord, a2.scene_coord)
            if (dmin is None or d < dmin) and d <= max_dist:
                dmin, amin1, amin2 = d, a1, a2
    report_closest(session, dmin, amin1, amin2, show, max_dist)
    return dmin, amin1, amin2

def register_command(session):
    from chimerax.core.commands import CmdDesc, register, FloatArg, BoolArg
    from chimerax.atomic import AtomsArg
    desc = CmdDesc(synopsis = 'find closest pair of atoms')
    register('humanize', desc, humanize, logger=session.logger)

register_command(session)
