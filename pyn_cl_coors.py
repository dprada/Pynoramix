import io_formats as io

class cl_traj():
    def __init__(self,file_input=None,frame=None,begin=None,end=None,increment=1,units=None,verbose=True):
        self.name=None
        self.type=None
        self.info=None
        self.io_file=None
        self.io_pos=None
        self.io_status=False
        self.io_vars=[]
        self.num_frames=0
        self.precision=None
        self.frame=[]

        if file_input!=None:
            self.name=file_input
            self.type=file_input.split('.')[1]
            self.io_status=False

        if frame!=None or begin!=None or end!=None:
            self.upload_frame(frame,begin,end,increment,units)

        if verbose:
            print '#',self.num_frames,'frames/models loaded.'

    def restart(self):
        self.io_status=None


    #def delete_frame(self,frame='ALL',begin=None,end=None,increment=1,units=None):
    # 
    #    if frame in ['all','All','ALL']:
    #        del(self.frame)
    #        self.frame=[]
    #        self.num_frames=0
    #        return
    # 
    #    elif type(frame) in [list,tuple]:
    #        for ii in frame:
    #            self.frame.__delitem__(ii)
    #            self.num_frames-=1
    #        return
    # 
    #    pass

    def reload_frame(self,old=0,new='next'):

        if new=='next':
            temp_frame,self.io_pos,ioerr=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos)
            self.frame[old]=temp_frame
            del(temp_frame)
        else:
            print '# Not supported yet.'
            return

    def upload_frame(self,frame='next',begin=None,end=None,increment=1,units=None):

        if not self.io_status:
            self.io_file,self.io_vars,self.io_pos,ioerr=getattr(io,'coor_'+self.type).open_traj(self.name)
            if ioerr: print '# Error opening the file'; return
            self.io_status=True


        if type(frame) in [int]: frame=[frame]

        if begin!=None or end!=None:
            if units==None:
                print "# Option units=None."
                print "# Choose units in ['ns','ps','md_steps','frames']"
                pass
            print "# Not supported yet."
            pass

        else:
            ### Uploading next frame
            if frame in ['next','Next','NEXT']:
                temp_frame,self.io_pos,ioerr=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos)
                if ioerr: return '# Error opening the file'
                self.frame.append(temp_frame)
                self.num_frames+=1
            ### Uploading all frames
            elif frame in ['all','All','ALL']:
                temp_frames,ioerr=getattr(io,'coor_'+self.type).read_all(self.io_file)
                if ioerr: return '# Error reading file'
                self.io_file,ioerr=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                if ioerr: return '# Error closing file'
                self.num_frames+=len(temp_frames)
                for ii in temp_frames:
                    self.frame.append(ii)
                self.precision=self.frame[0].precision
                self.io_status=False
            ### Uploading a list of frames
            elif type(frame) in [tuple,list]:
                for ii in frame:
                    temp_frame,self.io_pos,ioerr=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos,ii)
                    if ioerr: return '# Error reading file'
                    self.frame.append(temp_frame)
                    self.num_frames+=1

            else:
                print "# Options not supported yet."
                pass

        pass



        
