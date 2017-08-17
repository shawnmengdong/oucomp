function rock = setupRock(G)

rock = makeRock(G,100*milli*darcy, 0.13);
rock.cr = 4e-6/psia;
rock.pr = 3550*psia;

end